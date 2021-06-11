package dx.executor

import java.nio.file.Path
import dx.api.{DxApi, DxFile}

import scala.annotation.tailrec

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
    * @param listings mapping of source directory to Set of files/folders to
    *                 include in the upload; if a folder is included, all of its
    *                 files and subfolders are included regardless of whether
    *                 they appear in the Set individually
    * @return mapping of source path to (projectId, folder)
    */
  def uploadDirectories(dirs: Set[Path],
                        recursive: Boolean = true,
                        wait: Boolean = false,
                        listings: Map[Path, Set[Path]] = Map.empty): Map[Path, (String, String)]

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
    * @param listings mapping of source directory to Set of files/folders to
    *                 include in the upload; if a folder is included, all of its
    *                 files and subfolders are included regardless of whether
    *                 they appear in the Set individually
    * @return mapping of source path to (projectId, folder, files), where files
    *         is a Vector of all the uploaded files
    */
  def uploadDirectoriesWithDestination(
      dirs: Map[Path, String],
      recursive: Boolean = true,
      wait: Boolean = false,
      listings: Map[Path, Set[Path]] = Map.empty
  ): Map[Path, (String, String, Map[Path, DxFile])]
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

  private def includePath(path: Path, paths: Set[Path]): Boolean = {
    @tailrec
    def containsAncestor(child: Path): Boolean = {
      Option(child.getParent) match {
        case Some(parent) => paths.contains(parent) || containsAncestor(parent)
        case None         => false
      }
    }
    paths.contains(path) || (path.toFile.isFile && containsAncestor(path))
  }

  override def uploadDirectories(
      dirs: Set[Path],
      recursive: Boolean,
      wait: Boolean,
      listings: Map[Path, Set[Path]] = Map.empty
  ): Map[Path, (String, String)] = {
    dirs.map { path =>
      val filter = listings
        .get(path)
        .map(paths => (path: Path) => includePath(path, paths))
      val (project, folder, _) =
        dxApi.uploadDirectory(path, recursive = recursive, wait = wait, filter = filter)
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
  override def uploadDirectoriesWithDestination(
      dirs: Map[Path, String],
      recursive: Boolean,
      wait: Boolean,
      listings: Map[Path, Set[Path]] = Map.empty
  ): Map[Path, (String, String, Map[Path, DxFile])] = {
    dirs.map {
      case (path, destination) =>
        val filter = listings
          .get(path)
          .map(paths => (path: Path) => includePath(path, paths))
        val (project, folder, files) =
          dxApi.uploadDirectory(path, Some(destination), recursive = recursive, wait = wait, filter)
        path -> (project.getOrElse(dxApi.currentProject.id), folder, files)
    }
  }
}
