package dx.executor

import java.nio.file.Path
import dx.api.{DxApi, DxFile}

import scala.annotation.tailrec

/**
  * A file to upload.
  * @param source The source file.
  * @param destination Optional destination project and/or path; defaults to the
  *                    context project and root folder.
  * @param tags tags to add to the uploaded file.
  * @param properties properties to add to the uploaded file.
  */
case class FileUpload(source: Path,
                      destination: Option[String] = None,
                      tags: Set[String] = Set.empty,
                      properties: Map[String, String] = Map.empty)

/**
  * A directory to upload.
  * @param source The source directory
  * @param destination Optional destination project and/or folder; defaults to the
  *                    context project and root folder.
  * @param recursive Whether to search recursively in `source` for files/directories
  *                  to upload.
  * @param listing Set of files/folders to include in the upload; if a folder is
  *                included, all of its files and subfolders are included regardless
  *                of whether they appear in the Set individually
  */
case class DirectoryUpload(source: Path,
                           destination: Option[String] = None,
                           recursive: Boolean = true,
                           listing: Option[Set[Path]] = None)

trait FileUploader {

  /**
    * Uploads files.
    * @param files Set of `FileUpload` objects, one per file.
    * @param wait Whether to wait for uploads to complete before returning.
    * @return mapping of sorce Path to DxFile
    */
  def uploadFiles(files: Set[FileUpload], wait: Boolean = false): Map[Path, DxFile]

  /**
    * Uploads directories to the context project and folder.
    * @param dirs directories to upload
    * @param wait whether to wait for upload to complete
    * @return mapping of source path to (projectId, folder)
    */
  def uploadDirectories(dirs: Set[DirectoryUpload],
                        wait: Boolean = false): Map[Path, (String, String, Map[Path, DxFile])]
}

/**
  * Simple FileUploader that uploads one file at a time
  * using the API.
  * @param dxApi DxApi
  */
case class SerialFileUploader(dxApi: DxApi = DxApi.get) extends FileUploader {
  def uploadFiles(files: Set[FileUpload], wait: Boolean = false): Map[Path, DxFile] = {
    files.map {
      case FileUpload(path, dest, tags, properties) =>
        path -> dxApi.uploadFile(path, dest, wait = wait, tags = tags, properties = properties)
    }.toMap
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

  /**
    * Uploads directories to the context project and folder.
    *
    * @param dirs directories to upload
    * @param wait whether to wait for upload to complete
    * @return mapping of source path to (projectId, folder)
    */
  override def uploadDirectories(dirs: Set[DirectoryUpload],
                                 wait: Boolean): Map[Path, (String, String, Map[Path, DxFile])] = {
    dirs.map { dir =>
      val filter = dir.listing.map(paths => (path: Path) => includePath(path, paths))
      val (project, folder, files) =
        dxApi.uploadDirectory(dir.source,
                              dir.destination,
                              recursive = dir.recursive,
                              wait = wait,
                              filter)
      dir.source -> (project.getOrElse(dxApi.currentProjectId.get), folder, files)
    }.toMap
  }
}
