package dx.executor

import dx.api.{DxAnalysis, DxApi, DxFile, DxJob, FileUpload}
import dx.core.Constants
import dx.core.io.DxWorkerPaths
import dx.core.ir.{DxName, DxNameFactory, ParameterLinkSerializer, TypeSerde}
import dx.core.languages.wdl.{WdlBundle, WdlDxName, WdlUtils}
import dx.util.{CodecUtils, Logger, PosixPath}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsArray, JsNumber, JsString, JsValue}
import wdlTools.types.{TypedAbstractSyntax => TAT}

import java.nio.file.{Files, Path, Paths}

private object TestWorkerJobMeta {
  val InstanceType: String = "mem_ssd_unicorn"
}

private case class TestWorkerJobMeta(override val workerPaths: DxWorkerPaths,
                                     override val waitOnUpload: Boolean = false,
                                     override val dxNameFactory: DxNameFactory,
                                     override val dxApi: DxApi = DxApi.get,
                                     override val logger: Logger = Logger.get)
    extends WorkerJobMeta(workerPaths, waitOnUpload, dxNameFactory, dxApi, logger) {

  override val jobId: String = null

  override def codeFile: Path = {
    workerPaths.getRootDir().resolve(s"${jobId}.code.sh").asJavaPath
  }

}

class WorkerJobMetaTest extends AnyFlatSpec with Matchers {

  private def setup(): DxWorkerPaths = {
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()
    val workerPaths = DxWorkerPaths(PosixPath(jobRootDir.toString))
    workerPaths.createCleanDirs()
    workerPaths
  }

  // Note: if the file doesn't exist, this throws a null pointer exception
  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  private def parse(path: Path): WdlBundle = {
    val (doc, _) = WdlUtils.parseAndCheckSourceFile(path)
    WdlBundle.create(doc)
  }

  it should "fail at the internal python process and output subprocess error stack" in {
    val workerPaths = setup()
    val path = pathFromBasename("fixtures", "python_task_fail.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }

    val workerJobMeta = WorkerJobMeta(
        workerPaths = workerPaths,
        waitOnUpload = false,
        dxNameFactory = WdlDxName,
        dxApi = DxApi.createAndSet()
    )
    workerJobMeta.runJobScriptFunction(name = "test_script")

    workerJobMeta should be
    wf should be(None)
  }

}
