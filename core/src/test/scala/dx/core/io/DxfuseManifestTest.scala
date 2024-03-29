package dx.core.io

import java.nio.file.{Files, Path}
import dx.core.Assumptions.isLoggedIn
import dx.core.Tags.ApiTest
import dx.api.{DxApi, DxFile, DxProject}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import dx.util.{Logger, PosixPath}

class DxfuseManifestTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)
  private val dxApi: DxApi = DxApi()(Logger.Quiet)
  private val rootDir = Files.createTempDirectory("root")
  rootDir.toFile.deleteOnExit()
  private val dxPathConfig = DxWorkerPaths(PosixPath(rootDir.toString))
  private val ArchivedProj = "ArchivedStuff"
  private lazy val dxArchivedProj: DxProject = dxApi.resolveProject(ArchivedProj)

  it should "detect and provide legible error for archived files" taggedAs ApiTest in {
    val fileDir: Map[String, Path] = Map(
        s"dx://${ArchivedProj}:/Catch22.txt" -> dxPathConfig
          .getDxfuseMountDir()
          .resolve("inputs/A")
          .asJavaPath,
        s"dx://${ArchivedProj}:/LICENSE" -> dxPathConfig
          .getDxfuseMountDir()
          .resolve("inputs/B")
          .asJavaPath,
        s"dx://${ArchivedProj}:/README" -> dxPathConfig
          .getDxfuseMountDir()
          .resolve("inputs/C")
          .asJavaPath
    )

    // resolve the paths
    val (uris, resolvedFiles) =
      dxApi
        .resolveDataObjectBulk(fileDir.keys.toVector, dxArchivedProj)
        .map {
          case (dxUri, dxFile: DxFile) => (dxUri, dxFile)
          case other                   => throw new Exception(s"expected file, not ${other}")
        }
        .unzip

    // describe the files
    val describedFiles = dxApi.describeFilesBulk(resolvedFiles.toVector)
    val filesInManifest: Map[DxFile, Path] = uris
      .zip(describedFiles)
      .map {
        case (dxUri, dxFile: DxFile) => dxFile -> fileDir(dxUri)
      }
      .toMap

    // Creating a manifest should fail, because some of the files are archived
    assertThrows[Exception] {
      DxfuseManifestBuilder(dxPathConfig, dxApi).apply(filesInManifest, Map.empty)
    }
  }
}
