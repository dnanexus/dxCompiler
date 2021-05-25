package dx.executor

import java.nio.file.Path

import dx.api.{DxApi, DxFile}

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
