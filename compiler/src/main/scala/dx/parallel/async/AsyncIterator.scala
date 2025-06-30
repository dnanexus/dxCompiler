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

import java.util.concurrent.{ArrayBlockingQueue, BlockingQueue, LinkedBlockingDeque, TimeUnit}

object AsyncIterator {

  /** Builds an [[AsyncIterator]] and starts it.
    *
    * @param source the source iterator to consume
    * @param bufferSize the maximum number of elements to buffer before the blocking when consuming the source, or None
    *                   if unbounded.
    * @tparam T the type of object to consume
    */
  def apply[T](source: Iterator[T], bufferSize: Option[Int] = None): AsyncIterator[T] = {
    new AsyncIterator[T](source, bufferSize).start()
  }
}

/** An asynchronous wrapper for an [[Iterator]].  A separate thread will be created to consume the source iterator.
  * Will buffer up to [[bufferSize]] elements before the blocking when consuming the source.
  *
  * @param source the source iterator to consume
  * @param bufferSize the maximum number of elements to buffer before the blocking when consuming the source, or None
  *                   if unbounded.
  * @tparam T the type of object to consume
  */
class AsyncIterator[T](private val source: Iterator[T], bufferSize: Option[Int] = None)
    extends Iterator[T]
    with AsyncRunnable {
  bufferSize.foreach(size =>
    require(size > 0, s"bufferSize must be greater than zero when given, found $size")
  )

  private var buffer: Option[T] = None

  private val queue: BlockingQueue[T] = bufferSize match {
    case Some(size) => new ArrayBlockingQueue[T](size)
    case None       => new LinkedBlockingDeque[T]()
  }

  protected def execute(): Unit = this.source.foreach(this.queue.put)

  /** Returns true if there exists more elements, false otherwise */
  def hasNext: Boolean = {
    checkAndRaise()

    // Get the next item, or wait until the underlying thread is done and there are no more items in the queue
    while (buffer.isEmpty && !(this.done && this.queue.isEmpty)) {
      checkAndRaise() // check if the underlying thread raised an exception
      tryAndModifyInterruptedException("Interrupted waiting on taking from the queue.") {
        buffer = Option(this.queue.poll(50, TimeUnit.MILLISECONDS))
      }
    }

    // If there are no more elements in the source iterator then await the signal for final completion of the thread's
    // execution work. We must await the finished countdown latch first because a race condition may occur when the
    // source iterator raises an exception but we have not yet had a chance to save the exception message before
    // finishing the final call to `hasNext()`. Blocking and awaiting the final countdown latch signal guarantees
    // we will see the exception message if it is present. If it is present, then we get a final chance to re-raise it.
    // If we did not handle exceptions this way, then the iterator may be silently cut short and data truncated.
    // Issue: https://github.com/fulcrumgenomics/commons/pull/74
    if (buffer.isEmpty) {
      awaitDone()
      checkAndRaise()
    }
    buffer.nonEmpty
  }

  /** Gets the next item. */
  def next(): T = {
    checkAndRaise()
    this.buffer match {
      case None => throw new NoSuchElementException("Calling next() when hasNext() is false.")
      case Some(item) =>
        this.buffer = None
        item
    }
  }
}
