package dx.executor

import java.nio.file.Path

import dx.api.{DxApi, DxFile}

trait FileUploader {

  /**
    * Uploads files to the context project and folder.
    * @param files files to upload
    * @param wait whether to wait for upload to complete
    * @return
    */
  def uploadFiles(files: Set[Path], wait: Boolean = false): Map[Path, DxFile]

  /**
    * Uploads directories to the context project and folder.
    * @param dirs directories to upload
    * @param recursive whether to upload directories recursively
    * @param wait whether to wait for upload to complete
    * @return mapping of source path to (projectId, folder)
    */
  def uploadDirectories(dirs: Set[Path],
                        recursive: Boolean = true,
                        wait: Boolean = false): Map[Path, (String, String)]

  /**
    * Uploads files/directories to specific destinations.
    * @param files mapping of source path to destination path
    * @param wait whether to wait for upload to complete
    * @return
    */
  def uploadFilesWithDestination(files: Map[Path, String], wait: Boolean = false): Map[Path, DxFile]

  /**
    * Uploads directories to specific destinations.
    * @param dirs mapping of directory source path to destination path
    * @param recursive whether to upload directories recursively
    * @param wait whether to wait for upload to complete
    * @return mapping of source path to (projectId, folder)
    */
  def uploadDirectoriesWithDestination(dirs: Map[Path, String],
                                       recursive: Boolean = true,
                                       wait: Boolean = false): Map[Path, (String, String)]
}

/**
  * Simple FileUploader that uploads one file at a time
  * using the API.
  * @param dxApi DxApi
  */
case class SerialFileUploader(dxApi: DxApi = DxApi.get) extends FileUploader {
  def uploadFiles(files: Set[Path], wait: Boolean = false): Map[Path, DxFile] = {
    files.map(path => path -> dxApi.uploadFile(path, wait = wait)).toMap
  }

  override def uploadDirectories(dirs: Set[Path],
                                 recursive: Boolean,
                                 wait: Boolean): Map[Path, (String, String)] = {
    dirs.map { path =>
      val (project, folder, _) = dxApi.uploadDirectory(path, recursive = recursive, wait = wait)
      path -> (project.getOrElse(dxApi.currentProject.id), folder)
    }.toMap
  }

  def uploadFilesWithDestination(files: Map[Path, String],
                                 wait: Boolean = false): Map[Path, DxFile] = {
    files.map {
      case (path, destination) =>
        path -> dxApi.uploadFile(path, Some(destination), wait = wait)
    }
  }

  /**
    * Uploads directories to specific destinations.
    *
    * @param dirs      mapping of directory source path to destination path
    * @param recursive whether to upload directories recursively
    * @param wait      whether to wait for upload to complete
    * @return mapping of source path to (projectId, folder)
    */
  override def uploadDirectoriesWithDestination(dirs: Map[Path, String],
                                                recursive: Boolean,
                                                wait: Boolean): Map[Path, (String, String)] = {
    dirs.map {
      case (path, destination) =>
        val (project, folder, _) =
          dxApi.uploadDirectory(path, Some(destination), recursive = recursive, wait = wait)
        path -> (project.getOrElse(dxApi.currentProject.id), folder)
    }
  }
}
