// Create a manifest for download agent (https://github.com/dnanexus/dxda).
//
package dx.core.io

import java.nio.file.{Path, Paths}
import dx.api._
import dx.util.Logger
import spray.json._

case class DxdaManifest(value: JsObject)

case class DxdaManifestBuilder(dxApi: DxApi, logger: Logger = Logger.get) {
  /*
  Start with paths like this:
    "dx://dxCompiler_playground:/test_data/fileB",
    "dx://dxCompiler_playground:/test_data/fileC",

  and generate a manifest like this:

   val fileA = JsObject("id" -> JsString("file-FGqFGBQ0ffPPkYP19gBvFkZy"),
                         "name" -> JsString("fileA"),
                         "folder" -> JsString("/test_data"))

  JsObject(
      projId1 -> JsArray(fileA, ...),
      projId2 -> JsArray(...)
  )

   */

  // create a manifest entry for a single file
  private def createFileEntry(dxFile: DxFile, destination: Path): JsValue = {
    val destinationFile = destination.toFile
    val name = destinationFile.getName
    val folder = destinationFile.getParent
    val parts = dxFile
      .describe()
      .parts
      .map { parts =>
        Map("parts" -> JsObject(parts.map {
          case (idx, part) =>
            idx.toString -> JsObject(
                "size" -> JsNumber(part.size),
                "md5" -> JsString(part.md5)
            )
        }))
      }
      .getOrElse(Map.empty)
    JsObject(
        Map(
            "id" -> JsString(dxFile.id),
            "name" -> JsString(name),
            "folder" -> JsString(folder)
        ) ++ parts
    )
  }

  private def createFolderEntry(projectId: String,
                                folder: String,
                                destination: Path): Vector[JsValue] = {
    // dxda manifest doesn't support folders so we have to list the folder contents and add
    // all the files to the manifest
    val findDataObjects = DxFindDataObjects(dxApi)
    val project = dxApi.project(projectId)
    val result = findDataObjects
      .apply(
          Some(project),
          Some(folder),
          recurse = true,
          classRestriction = Some("file"),
          state = Some(DxState.Closed),
          extraFields = Set(Field.Parts)
      )
      .keys
      .toVector
    // logger.trace(s"in ${projectId}:${folder} found: ${result}")
    val folderPath = Paths.get(folder)
    result.map {
      case dxFile: DxFile =>
        val fileRelFolder = folderPath.relativize(Paths.get(dxFile.getFolder))
        val fileDest = destination.resolve(fileRelFolder).resolve(dxFile.getName)
        createFileEntry(dxFile, fileDest)
      case other =>
        throw new Exception(s"not a file: ${other}")
    }
  }

  /**
    *
    * @param fileToLocalMapping mapping of
    * @return
    */
  def apply(fileToLocalMapping: Map[DxFile, Path],
            folderToLocalMapping: Map[(String, String), Path]): Option[DxdaManifest] = {
    if (fileToLocalMapping.isEmpty && folderToLocalMapping.isEmpty) {
      return None
    }

    val filesByContainer: Map[DxProject, Vector[DxFile]] =
      fileToLocalMapping.keys.toVector.groupBy { file =>
        // make sure file is pre-described
        assert(file.hasCachedDesc)
        // make sure file is in the live state - archived files cannot be accessed.
        if (file.describe().archivalState != DxArchivalState.Live) {
          throw new Exception(s"file ${file} is not live")
        }
        // create a sub-map per container
        dxApi.project(file.describe().project)
      }

    val foldersByContainer: Map[DxProject, Vector[String]] =
      folderToLocalMapping.keys.toVector.groupBy(_._1).map {
        case (projectId, folders) => dxApi.project(projectId) -> folders.map(_._2)
      }

    val manifest: Map[String, JsValue] =
      (filesByContainer.keySet ++ foldersByContainer.keySet).map { dxContainer =>
        val containerFiles = filesByContainer.getOrElse(dxContainer, Vector.empty)
        val containerFilesToLocalPath: Vector[JsValue] =
          containerFiles.map { dxFile =>
            createFileEntry(dxFile, fileToLocalMapping(dxFile))
          }
        val containerFolders = foldersByContainer.getOrElse(dxContainer, Vector.empty)
        val containerFolderFilesToLocalPath: Vector[JsValue] =
          containerFolders.flatMap { folder =>
            createFolderEntry(dxContainer.id,
                              folder,
                              folderToLocalMapping((dxContainer.id, folder)))
          }
        dxContainer.id -> JsArray(containerFilesToLocalPath ++ containerFolderFilesToLocalPath)
      }.toMap

    Some(DxdaManifest(JsObject(manifest)))
  }
}
