package dx.executor

import java.nio.file.Path
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
    * @return
    */
  def upload(files: Set[FileUpload]): Map[Path, DxFile]
}

object FileUploader {
  //      /**
  //        * Times the execution of a block of code.
  //        * TODO: move this to dxCommon SysUtils
  //        * @param block the block to execute
  //        * @tparam R the return value type
  //        * @return (return value, time in seconds)
  //        */
  def time[R](timeUnit: TimeUnit = TimeUnit.SECONDS)(block: => R): (R, Long) = {
    val t0 = System.nanoTime()
    val result = block
    val t1 = System.nanoTime()
    (result, timeUnit.convert(t1 - t0, TimeUnit.NANOSECONDS))
  }
}

/**
  * Simple FileUploader that uploads one file at a time
  * using the API.
  * @param dxApi DxApi
  */
case class SerialFileUploader(waitOnUpload: Boolean = false,
                              dxApi: DxApi = DxApi.get,
                              logger: Logger = Logger.get)
    extends FileUploader {
  def upload(files: Set[FileUpload]): Map[Path, DxFile] = {
    val (results, duration) = FileUploader.time() {
      files.map {
        case FileUpload(path, dest, tags, properties) =>
          path -> dxApi.uploadFile(path,
                                   dest,
                                   wait = waitOnUpload,
                                   tags = tags,
                                   properties = properties)
      }.toMap
    }
    logger.trace(s"Uploaded ${results.size} files in ${duration} seconds")
    results
  }
}

case class ParallelFileUploader(waitOnUpload: Boolean = false,
                                maxConcurrent: Int = Runtime.getRuntime.availableProcessors(),
                                dxApi: DxApi = DxApi.get,
                                logger: Logger = Logger.get)
    extends FileUploader {

  private case class UploadCallable(upload: FileUpload) extends Callable[(Path, DxFile)] {
    override def call(): (Path, DxFile) = {
      upload.source -> dxApi.uploadFile(upload.source,
                                        upload.destination,
                                        wait = waitOnUpload,
                                        tags = upload.tags,
                                        properties = upload.properties)
    }
  }

  def upload(files: Set[FileUpload]): Map[Path, DxFile] = {
    val callables = files.map(UploadCallable).toVector
    if (callables.size == 1) {
      Map(callables.head.call())
    } else {
      val executor = Executors.newFixedThreadPool(Math.min(files.size, maxConcurrent))

      def shutdown(now: Boolean = false): Unit = {
        try {
          if (now) {
            executor.shutdownNow()
          } else {
            executor.shutdown()
          }
        } catch {
          case se: SecurityException =>
            logger.warning(
                "Unexpected security exception shutting down upload thread pool executor",
                exception = Some(se)
            )
        }
      }

      // submit all upload jobs
      val futures =
        try {
          callables.map { callable: Callable[(Path, DxFile)] =>
            executor.submit(callable)
          }
        } catch {
          case re: RejectedExecutionException =>
            shutdown(now = true)
            throw new Exception("Unexpected rejection of file upload task", re)
        }

      // try to shut down the threadpool gracefully - this
      // will let the submitted jobs finish
      shutdown()

      // wait for all jobs to complete - throw exception if any of
      // the uploads fail
      val (results, duration) = FileUploader.time() {
        futures.zipWithIndex.map {
          case (f, index) =>
            try {
              f.get
            } catch {
              case t: Throwable =>
                shutdown(now = true)
                throw new Exception(s"Error uploading ${callables(index).upload.source}", t)
            }
        }
      }

      logger.trace(s"Uploaded ${results.size} files in ${duration} seconds")

      results.toMap
    }
  }
}
