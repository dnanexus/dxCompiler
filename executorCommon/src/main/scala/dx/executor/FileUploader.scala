package dx.executor

import java.nio.file.{Files, Path}
import java.util.concurrent.{Callable, Executors, RejectedExecutionException, TimeUnit}

import dx.api.{DxApi, DxFile}
import dx.util.Logger

case class FileUpload(source: Path,
                      destination: Option[String] = None,
                      tags: Set[String] = Set.empty,
                      properties: Map[String, String] = Map.empty)

trait FileUploader {

  /**
    * Upload files to the context project and folder.
    * @param files files to upload
    * @param wait whether to wait for upload to complete
    * @return
    */
  def upload(files: Set[FileUpload], wait: Boolean = false): Map[Path, DxFile]
}

/**
  * Simple FileUploader that uploads one file at a time
  * using the API.
  * @param dxApi DxApi
  */
case class SerialFileUploader(dxApi: DxApi = DxApi.get) extends FileUploader {
  def upload(files: Set[FileUpload], wait: Boolean = false): Map[Path, DxFile] = {
    files.map {
      case FileUpload(path, dest, tags, properties) =>
        path -> dxApi.uploadFile(path, dest, wait = wait, tags = tags, properties = properties)
    }.toMap
  }
}

case class ParallelFileUploader(maxConcurrent: Int = Runtime.getRuntime.availableProcessors(),
                                dxApi: DxApi = DxApi.get,
                                logger: Logger = Logger.get)
    extends FileUploader {
  private case class UploadCallable(upload: FileUpload, waitOnUpload: Boolean)
      extends Callable[(Path, DxFile)] {
    override def call(): (Path, DxFile) = {
      upload.source -> dxApi.uploadFile(upload.source,
                                        upload.destination,
                                        wait = waitOnUpload,
                                        tags = upload.tags,
                                        properties = upload.properties)
    }
  }

  def upload(files: Set[FileUpload], wait: Boolean = false): Map[Path, DxFile] = {
    val callables = files.map(UploadCallable(_, wait)).toVector
    if (callables.size == 1) {
      Map(callables.head.call())
    } else {
      val executor = Executors.newFixedThreadPool(maxConcurrent)

//      /**
//        * Times the execution of a block of code.
//        * TODO: move this to dxCommon SysUtils
//        * @param block the block to execute
//        * @tparam R the return value type
//        * @return (return value, time in seconds)
//        */
      def time[R](timeUnit: TimeUnit = TimeUnit.SECONDS)(block: => R): (R, Long) = {
        val t0 = System.nanoTime()
        val result = block // call-by-name
        val t1 = System.nanoTime()
        (result, timeUnit.convert(t1 - t0, TimeUnit.NANOSECONDS))
      }

      def shutdown(force: Boolean = false): (Boolean, Long) = {
        val timeoutSecs = if (force) {
          10L
        } else {
          // Upload could take a long time, especially if we're waiting for them
          // to complete. If so, we base the timeout on the file size.
          if (wait) {
            Math.round(files.map(f => Files.size(f.source)).sum / 10_000_000L)
          } else {
            files.size * 10L
          }
        }
        try {
          if (force) {
            executor.shutdownNow()
          } else {
            executor.shutdown()
          }
          time() { executor.awaitTermination(timeoutSecs, TimeUnit.SECONDS) }
        } catch {
          case se: SecurityException =>
            logger.warning(
                "Unexpected security exception shutting down upload thread pool executor",
                exception = Some(se)
            )
            (false, timeoutSecs)
          case ie: InterruptedException =>
            logger.warning("Interrupted while awaiting upload thread pool termination",
                           exception = Some(ie))
            (false, timeoutSecs)
        }
      }

      val futures =
        try {
          callables.map(executor.submit)
        } catch {
          case re: RejectedExecutionException =>
            shutdown(force = true)
            throw new Exception("Unexpected rejection of file upload task", re)
        }

      // wait for all jobs to complete
      val (terminated, duration) = shutdown()
      if (terminated) {
        logger.trace(s"Upload of ${files.size} files completed in ${duration} seconds")
      } else {
        throw new Exception(
            s"Uploading ${files.size} files did not complete within ${duration} seconds"
        )
      }

      val results = futures.zipWithIndex.map {
        case (f, index) =>
          try {
            Some(f.get)
          } catch {
            case t: Throwable =>
              val callable = callables(index)
              logger.error(s"Error uploading ${callable.upload.source}", Some(t))
              None
          }
      }

      val (successes, failures) = results.partition(_.isDefined)
      if (failures.nonEmpty) {
        throw new Exception(s"Upload failed for ${failures.size} files")
      }
      successes.flatten.toMap
    }
  }
}
