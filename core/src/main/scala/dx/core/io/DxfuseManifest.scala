// Create a manifest for dxfuse (https://github.com/dnanexus/dxfuse).
//

package dx.core.io

import java.nio.file.Path
import dx.api.{DxApi, DxArchivalState, DxFile}
import dx.util.{Logger, PosixPath}
import spray.json._

case class DxfuseManifest(value: JsValue)

case class DxfuseManifestBuilder(workerPaths: DxWorkerPaths,
                                 dxApi: DxApi,
                                 logger: Logger = Logger.get) {
  def apply(
      fileToLocalMapping: Map[DxFile, Path] = Map.empty,
      folderToLocalMapping: Map[(String, String), Path] = Map.empty,
      folderListings: Map[(String, String), Set[PosixPath]] = Map.empty
  ): Option[DxfuseManifest] = {
    if (fileToLocalMapping.isEmpty) {
      return None
    }

    val files = fileToLocalMapping.map {
      case (dxFile, path) =>
        // we expect that the files will have already been bulk described
        assert(dxFile.hasCachedDesc)
        // check that the files are not archived
        if (dxFile.describe().archivalState != DxArchivalState.Live) {
          throw new Exception(s"file ${dxFile} is not live")
        }

        val parentDir = path.getParent
        // remove the mountpoint from the directory. We need
        // paths that are relative to the mount point.
        val mountDir = workerPaths.getDxfuseMountDir().asJavaPath
        assert(parentDir.startsWith(mountDir))
        val relParentDir = s"/${mountDir.relativize(parentDir)}"

        val desc = dxFile.describe()
        JsObject(
            "file_id" -> JsString(dxFile.id),
            "parent" -> JsString(relParentDir),
            "proj_id" -> JsString(desc.project),
            "fname" -> JsString(path.getFileName.toString),
            "size" -> JsNumber(desc.size),
            "ctime" -> JsNumber(desc.created),
            "mtime" -> JsNumber(desc.modified)
        )
    }.toVector

    val folders = folderToLocalMapping.map {
      case ((project, folder), path) =>
        JsObject(
            "proj_id" -> JsString(project),
            "folder" -> JsString(folder),
            "dirname" -> JsString(path.toString)
        )
    }.toVector

    Some(
        DxfuseManifest(
            JsObject("files" -> JsArray(files), "directories" -> JsArray(folders))
        )
    )
  }
}
