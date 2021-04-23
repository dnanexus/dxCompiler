package dx.executor

import java.nio.file.Path

import dx.api.{DxApi, DxFile}

trait FileUploader {

  /**
    * Upload files to the context project and folder.
    * @param files files to upload
    * @param wait whether to wait for upload to complete
    * @return
    */
  def upload(files: Set[Path], wait: Boolean = false): Map[Path, DxFile]

  /**
    * Upload files to specific destinations.
    * @param files mapping of source file to destination path
    * @param wait whether to wait for upload to complete
    * @return
    */
  def uploadWithDestination(files: Map[Path, String], wait: Boolean = false): Map[Path, DxFile]
}

/**
  * Simple FileUploader that uploads one file at a time
  * using the API.
  * @param dxApi DxApi
  */
case class SerialFileUploader(dxApi: DxApi = DxApi.get) extends FileUploader {
  def upload(files: Set[Path], wait: Boolean = false): Map[Path, DxFile] = {
    files.map(path => path -> dxApi.uploadFile(path, wait = wait)).toMap
  }

  def uploadWithDestination(files: Map[Path, String], wait: Boolean = false): Map[Path, DxFile] = {
    files.map {
      case (path, destination) =>
        (path, dxApi.uploadFile(path, Some(destination), wait = wait))
    }
  }
}
