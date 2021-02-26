package dx.executor

import java.nio.file.Path

import dx.api.{DxApi, DxFile}

trait FileUploader {

  /**
    * Upload files to the context project and folder.
    * @param files files to upload
    * @return
    */
  def upload(files: Set[Path]): Map[Path, DxFile]

  /**
    * Upload files to specific destinations.
    * @param files mapping of source file to destination path
    * @return
    */
  def upload(files: Map[Path, String]): Map[Path, DxFile]
}

/**
  * Simple FileUploader that uploads one file at a time
  * using the API.
  * @param dxApi DxApi
  */
case class SerialFileUploader(dxApi: DxApi = DxApi.get) extends FileUploader {
  def upload(files: Set[Path]): Map[Path, DxFile] = {
    files.map(path => path -> dxApi.uploadFile(path)).toMap
  }

  def upload(files: Map[Path, String]): Map[Path, DxFile] = {
    files.map {
      case (path, destination) =>
        (path, dxApi.uploadFile(path, Some(destination)))
    }
  }
}
