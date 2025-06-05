/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 Fulcrum Genomics LLC
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

package dx.parallel

import dx.parallel.collection.ParIterator

import scala.collection.parallel.immutable
import scala.collection.parallel.{ForkJoinTaskSupport, ParIterable, TaskSupport}
import java.util.concurrent.ForkJoinPool
import scala.language.implicitConversions

trait ParallelDef {

  /**
    * Implicit that provides additional methods to any collection that is Parallelizable.
    * Introduces [[parWith()]] methods that create parallel versions of the collection
    * with various configuration options.
    *
    * @param sequential a sequential (i.e. non-parallel) iterable
    * @param f a function that generates a parallel collection from the sequential collection
    * @tparam A the type of the elements in the collection
    * @tparam S the type of the non-parallel collection
    * @tparam P the type of the parallel collection
    *
    */
  class ParSupport[A, S <: Iterable[A], P <: ParIterable[A]](private val sequential: S, f: S => P) {

    /** Creates a parallel collection with the provided TaskSupport. */
    def parWith(taskSupport: TaskSupport): P = {
      val par = f(sequential)
      par.tasksupport = taskSupport
      par
    }

    /** Creates a parallel collection with the provided ForkJoinPool. */
    def parWith(pool: ForkJoinPool): P = parWith(taskSupport = new ForkJoinTaskSupport(pool))

    /** Creates a parallel collection with the desired level of parallelism and FIFO semantics. */
    def parWith(parallelism: Int, fifo: Boolean = true): P = {
      parWith(
          new ForkJoinPool(parallelism, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, fifo)
      )
    }
  }

  /** Implicit class that allows generation of ParIterators from Iterators with convenience methods. */
  implicit class ParIteratorSupport[A](iterator: Iterator[A]) {

    /**
      * Creates a [[ParIterator]]; see documentation for that class for detailed usage.
      *
      * @param threads the number of threads to use in parallel transform operations on the iterator
      * @param chunkSize the size of chunks to collect and perform parallel operations on
      * @param chunkBuffer if > 0 use an [[dx.parallel.async.AsyncIterator]] to accumulate/cache
      *                    `chunkBuffer` incoming chunks ready for parallel processing
      */
    def parWith(threads: Int,
                chunkSize: Int = ParIterator.DefaultChunkSize,
                chunkBuffer: Int = ParIterator.DefaultChunkBuffer): ParIterator[A] = {
      ParIterator(this.iterator,
                  threads = threads,
                  chunkSize = chunkSize,
                  chunkBuffer = chunkBuffer)
    }
  }

  /** Implicit that generates a ParSupport from a Seq. */
  implicit def seqToParSupport[A](
      seq: Seq[A]
  ): ParSupport[A, _ <: Seq[A], _ <: immutable.ParSeq[A]] = seq match {
    case v: Vector[A] => new ParSupport(v, (x: Vector[A]) => new immutable.ParVector[A](x))
    case s            => new ParSupport(s, (x: Seq[A]) => new immutable.ParVector[A](x.toVector))
  }
}

/** A singleton object providing access to all the functionality of ParallelDef. */
object ParallelDef extends ParallelDef
