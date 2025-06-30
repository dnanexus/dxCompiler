/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package dx.parallel.async

import java.util.concurrent.CountDownLatch
import java.util.concurrent.atomic.{AtomicInteger, AtomicReference}

import scala.concurrent.{Awaitable, ExecutionContext, Future}

object AsyncRunnable {

  /** A counter used for thread naming in [[AsyncRunnable]] and sub-classes */
  private val threadsCreated = new AtomicInteger(0)

  /** Generates a name for thread that will run the given runnable. */
  private[AsyncRunnable] def nextName(asyncRunnable: AsyncRunnable): String = {
    asyncRunnable.getClass.getSimpleName + this.threadsCreated.incrementAndGet
  }
}

/** A trait that can be mixed in to help manage asynchronous execution while being [[Runnable]].
  *
  * The `execute()` method should contain the work to be performed asynchronously.
  *
  * The `thread()` method can be used to return a [[Thread]] that wraps this [[Runnable]].
  *
  * The `start()` method can be used to create and start daemon [[Thread]] that wraps this [[Runnable]].
  *
  * The `checkAndRaise()` method should be used by implementing classes to check if the
  * `execute()`  method has encountered an exception, and if so, it will rethrow the exception as an
  * [[Error]] or [[RuntimeException]].
  *
  * The `tryAndModifyInterruptedException()` method should wrap blocking code in implementing classes
  * that are not part of the `execute()` method.  This method will throw a [[RuntimeException]] if interrupted while
  * waiting.
  */
trait AsyncRunnable extends Runnable {

  /** A reference to an exception thrown during an asynchronous operation. */
  private[AsyncRunnable] val _throwable: AtomicReference[Throwable] = new AtomicReference(null)

  /** This signals that the run method has started. */
  private val startedLatch: CountDownLatch = new CountDownLatch(1)

  /** This signals that the run method has completed, even if an exception occurred. This is needed so the close method
    * can be sure that the run method is no longer executing. */
  private val doneLatch: CountDownLatch = new CountDownLatch(1)

  final def run(): Unit = {
    startedLatch.countDown()
    try {
      this.execute()
    } catch {
      case thr: Throwable =>
        this._throwable.compareAndSet(null, thr)
        this.uponException()
    } finally {
      this.uponFinally()
      this.doneLatch.countDown()
    }
  }

  /** Returns true if the [[run()]] method has started, false otherwise. */
  def started: Boolean = this.startedLatch.getCount == 0

  /** Returns true if the [[run()]] method has completed, false otherwise. */
  def done: Boolean = this.doneLatch.getCount == 0

  /** Returns a throwable if an exception occurred in the [[run()]] method, None otherwise. */
  def throwable: Option[Throwable] = Option(this._throwable.get())

  /** Waits for the [[run()]] method to start. */
  def awaitStart(): Unit = this.startedLatch.await()

  /** Returns an [[Awaitable]] that completes when the [[run()]] method has started.  Returns the throwable if an
    * exception was encountered, [[None]] otherwise. */
  def uponStart()(implicit ec: ExecutionContext): Awaitable[Unit] = Future { this.awaitStart() }

  /** Waits for the [[run()]] method to complete. */
  def awaitDone(): Unit = this.doneLatch.await()

  /** Returns an [[Awaitable]] that completes when the [[run()]] method has completed.  Returns the throwable if an
    * exception was encountered, [[None]] otherwise. */
  def uponDone()(implicit ec: ExecutionContext): Awaitable[Option[Throwable]] = Future {
    this.awaitDone(); Option(this._throwable.get())
  }

  /** The name of this runnable.  This is used as the name of the thread in [[thread()]] as well.  The name is created
    * based on the class name and the number of [[AsyncRunnable]]s already created.
    * */
  val name: String = AsyncRunnable.nextName(this)

  /** Creates a new thread wrapping this runnable; the thread is not started. */
  final def thread(name: Option[String] = None): Thread = {
    require(!this.started, "Already started.")
    new Thread(this, name.getOrElse(this.name))
  }

  /** Starts this [[Runnable]] in a daemon thread.
    *
    * @param name optionally the name of the thread, otherwise a name is created based on the class name and the number
    *             of threads already created.
    */
  final def start(name: Option[String] = None, daemon: Boolean = true): this.type = {
    require(!this.started, "Already started.")
    val t = thread(name = name)
    t.setDaemon(daemon)
    t.start()
    this.startedLatch.countDown()
    this
  }

  /** The method that does the asynchronous work. */
  protected def execute(): Unit

  /** The method to execute if an exception occurs in the asynchronous thread.  This should not block. */
  protected def uponException(): Unit = ()

  /** The method to execute upon successfully execution of the run method or an exception occurs.  This should not block. */
  protected def uponFinally(): Unit = ()

  /** Checks to see if an exception has been raised by an asynchronous thread and if so rethrows it.  Use this method
    * before code that assumes the threads have not encountered an exception.
    */
  protected final def checkAndRaise(): Unit = {
    _throwable.getAndSet(null) match {
      case null           => ()
      case thr: Throwable => throw thr
    }
  }

  /** Executes the given block of code.  If an [[InterruptedException]] is thrown, throws a [[RuntimeException]]
    * with the given message.  Use this method for blocking code that when interrupted should not be recoverable.
    *
    * @param message the message to use if an [[InterruptedException]] is thrown by the block of code
    * @param f the block of code to execute
    */
  protected def tryAndModifyInterruptedException[T](message: String)(f: => T): T = {
    try {
      f
    } catch {
      case ie: InterruptedException => throw new RuntimeException(message, ie)
    }
  }
}
