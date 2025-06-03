/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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
 */

package dx.parallel.collection

import dx.parallel.ParallelDef._
import dx.parallel.async.AsyncIterator

import java.util.concurrent.ForkJoinPool
import scala.collection.AbstractIterator
import scala.collection.parallel.{ForkJoinTaskSupport, TaskSupport}

/**
  * Methods to help with manufacturing [[ParIterator]]s.
  */
object ParIterator {

  /** The default chunk size used by [[ParIterator]]s. */
  val DefaultChunkSize: Int = 2048

  /** The default incoming chunk buffering strategy for [[ParIterator]]s. */
  val DefaultChunkBuffer: Int = 2

  /**
    * Constructs a [[ParIterator]] from the underlying iterator that will perform operations in parallel on chunks
    * of `chunkSize` using `threads` threads.
    *
    * @param iter the underlying iterator to parallelize over
    * @param threads the number of threads to use for parallel computations like `map()`
    * @param chunkSize the number of elements to chunk together for parallel processing
    * @param chunkBuffer if > 0 use an  [[AsyncIterator]] to pre-buffer `chunkBuffer` chunks from the input
    *                    `iter` ready for parallel computation
    */
  def apply[A](iter: Iterator[A],
               threads: Int,
               chunkSize: Int = DefaultChunkSize,
               chunkBuffer: Int = DefaultChunkBuffer): ParIterator[A] = {
    val pool =
      new ForkJoinPool(threads, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true)
    val support = new ForkJoinTaskSupport(pool)
    apply(iter, support, chunkSize, chunkBuffer)
  }

  /**
    * Constructs a ParIterator from the underlying iterator, that will perform operations in parallel on chunks
    * of `chunkSize` using the provided `TaskSupport`.
    *
    * @param iter the underlying iterator to parallelize over
    * @param support the TaskSupport to use for parallel computations like `map()`
    * @param chunkSize the number of elements to chunk together for parallel processing
    * @param chunkBuffer if > 0 use an  [[AsyncIterator]] to pre-buffer `chunkBuffer` chunks from the input
    *                    `iter` ready for parallel computation
    */
  def apply[A](iter: Iterator[A],
               support: TaskSupport,
               chunkSize: Int,
               chunkBuffer: Int): ParIterator[A] = {
    val grouped = {
      val g = iter.grouped(chunkSize)
      if (chunkBuffer <= 0) g else AsyncIterator(g, Some(chunkBuffer))
    }

    new ParIterator[A](grouped, chunkSize, support)
  }
}

/**
  * A parallel iterator that operates by grouping or chunking the underlying iterator and then performing
  * operations on the chunks in parallel.  Optionally supports wrapping the iterator in an [[AsyncIterator]]
  * to draw items through the parallel processing and make them available for consumption while the next chunk(s)
  * are processed.
  *
  * For example:
  *
  * {{{
  * val iter = Io.readLines(reallyBigFile)
  * ParIterator(iter, chunkSize=4096, threads=8)
  *   .filter(_.contains("foo"))
  *   .map(_.toUpperCase)
  *   .map(_.reverse)
  *   .toAsync(cache=4096 * 8)
  *   .foreach(println)
  * }}}
  */
class ParIterator[A] private (private val iter: Iterator[Seq[A]],
                              val chunkSize: Int,
                              private val taskSupport: TaskSupport)
    extends AbstractIterator[A] {
  private var currentChunk: Iterator[A] = if (iter.hasNext) iter.next().iterator else Iterator.empty

  /** True if there are more items available, false otherwise. */
  override def hasNext: Boolean = {
    while (currentChunk.isEmpty && iter.hasNext) currentChunk = iter.next().iterator
    currentChunk.nonEmpty
  }

  /** Retrieve the next item.  Will throw an exception if there are no more items. */
  override def next(): A = {
    if (!hasNext) throw new NoSuchElementException
    currentChunk.next()
  }

  /** Returns an iterator over all remaining chunks including the current chunk. */
  private def chunks: Iterator[Seq[A]] = Iterator(currentChunk.to(LazyList)) ++ iter

  /** Wraps an async iterator around this parallel iterator to draw items through the iterator. This should only
    * be invoked after any other transforming operations such as map/filter/etc.
    *
    * @param cache how many elements to accumulate in the async iterator's cache - for best memory usage this
    *              should be a multiple of `chunkSize`.
    */
  def toAsync(cache: Int): AsyncIterator[A] = AsyncIterator(this, Some(cache))

  /** Wraps an async iterator around this parallel iterator to draw items through the iterator. This should only
    * be invoked after any other transforming operations such as map/filter/etc.
    */
  def toAsync: AsyncIterator[A] = toAsync(4 * chunkSize)
  // ^ Separate method so users can call .toAsync without () which is not doable with just a default above

  //////////////////////////////////////////////////////////////////////////////
  // Everything below here is just overrides of transform methods on Iterator
  //////////////////////////////////////////////////////////////////////////////

  override def map[B](f: A => B): ParIterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).map(f).seq), chunkSize, taskSupport)

  override def flatMap[B](f: A => IterableOnce[B]): ParIterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).flatMap(f).seq),
                       chunkSize,
                       taskSupport)

  override def filter(p: A => Boolean): ParIterator[A] =
    new ParIterator[A](chunks.map(c => c.parWith(taskSupport).filter(p).seq),
                       chunkSize,
                       taskSupport)

  override def filterNot(p: A => Boolean): ParIterator[A] =
    new ParIterator[A](chunks.map(c => c.parWith(taskSupport).filterNot(p).seq),
                       chunkSize,
                       taskSupport)

  override def withFilter(p: A => Boolean): ParIterator[A] =
    new ParIterator[A](chunks.map(c => c.parWith(taskSupport).withFilter(p).seq),
                       chunkSize,
                       taskSupport)

  override def collect[B](pf: PartialFunction[A, B]): ParIterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).collect(pf).seq),
                       chunkSize,
                       taskSupport)

  override def flatten[B](implicit ev: A => IterableOnce[B]): Iterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).flatten.seq), chunkSize, taskSupport)
}
