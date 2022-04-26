package dx.translator.cwl

import dx.core.ir.{Application, Callable, InstanceTypeSelection, Workflow}
import dx.translator.DefaultReorgSettings
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import dx.core.languages.cwl.CwlBundle
import dx.core.io.DxWorkerPaths
import dx.core.languages.Language
import dx.core.languages.cwl.DxHintSchema
import dx.cwl.{Process, Parser, ParserResult}
import dx.util.PosixPath
import java.nio.file.{Files, Path, Paths}

class ProcessTranslatorTest extends AnyFlatSpec with Matchers {

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  private def getProcessDocContent(processName: String): Map[String, String] = {
    val cwlFile: Path = pathFromBasename("cwl", s"${processName}.cwl.json")
    // val inputs = getInputs(cwlName)
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()
    val workerPaths = DxWorkerPaths(PosixPath(jobRootDir.toString))
    workerPaths.createCleanDirs()

    val parser = Parser.create(hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClassFromFile(cwlFile) match {
      case (version, _) if Language.parse(version) == Language.CwlV1_2 => ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL document of a supported version
               |${cwlFile}""".stripMargin
        )
    }
    val (process, schemas) = parser.parseFile(cwlFile) match {
      case ParserResult(Some(process: Process), _, _, schemas) => (process, schemas)
      case other =>
        throw new Exception(s"expected CWL document to contain a CommandLineTool, not ${other}")
    }
    val cwlBundle: CwlBundle = CwlBundle.create(process)
    val callableTranslator = ProcessTranslator(
        cwlBundle,
        cwlSchemas = schemas,
        locked = false,
        defaultRuntimeAttrs = Map.empty,
        reorgAttrs = DefaultReorgSettings(false),
        perWorkflowAttrs = Map.empty,
        defaultScatterChunkSize = 1,
        useManifests = false,
        instanceTypeSelection = InstanceTypeSelection.Static
    )

    val depOrder: Vector[Process] = cwlBundle.sortByDependencies

    // translate processes
    val (allCallables, sortedCallables) =
      depOrder.foldLeft((Map.empty[String, Callable], Vector.empty[Callable])) {
        case ((allCallables, sortedCallables), callable) =>
          val isPrimary = callable.name == cwlBundle.primaryProcess.name
          val translatedCallables =
            callableTranslator.translateProcess(callable, allCallables, isPrimary = isPrimary)
          (
              allCallables ++ translatedCallables.map(c => c.name -> c).toMap,
              sortedCallables ++ translatedCallables
          )
      }

    sortedCallables.map {
      case Application(name, _, _, _, _, _, document, _, _, _, _, _) =>
        name -> document.getDocContents
      case Workflow(name, _, _, _, document, _, _, _, _, _, _) =>
        name -> document.getDocContents
    }.toMap
  }

  "ProcessTranslator" should "render same cwl code for every unchanged process and different if there are changes" in {
    val standAloneProcessDoc1 = getProcessDocContent("count-lines1-wf")
    val standAloneProcessDoc2 = getProcessDocContent("count-lines1-wf-dup")
    standAloneProcessDoc1("parseInt-tool.cwl") should equal(
        standAloneProcessDoc2("parseInt-tool.cwl")
    )

    // task wc-tool.cwl changed
    standAloneProcessDoc1("wc-tool.cwl") should not equal standAloneProcessDoc2("wc-tool.cwl")
    standAloneProcessDoc1("count-lines1-wf_common") should not equal standAloneProcessDoc2(
        "count-lines1-wf_common"
    )
    standAloneProcessDoc1("count-lines1-wf_outputs") should not equal standAloneProcessDoc2(
        "count-lines1-wf_outputs"
    )
    standAloneProcessDoc1("count-lines1-wf") should not equal standAloneProcessDoc2(
        "count-lines1-wf"
    )

  }

  "ProcessTranslator" should "render different cwl code for changed common/outputs" in {

    val standAloneProcessDoc1 = getProcessDocContent("scatter-valuefrom-wf2")
    val standAloneProcessDoc2 = getProcessDocContent("scatter-valuefrom-wf2-dup")
    standAloneProcessDoc1("step1command") should equal(standAloneProcessDoc2("step1command"))

    // top-level wf input/output changed
    standAloneProcessDoc1("scatter-valuefrom-wf2_frag_stage-0") should not equal (
        standAloneProcessDoc2("scatter-valuefrom-wf2_frag_stage-0")
    )
    standAloneProcessDoc1("scatter-valuefrom-wf2_common") should not equal standAloneProcessDoc2(
        "scatter-valuefrom-wf2_common"
    )
    standAloneProcessDoc1("scatter-valuefrom-wf2_outputs") should not equal standAloneProcessDoc2(
        "scatter-valuefrom-wf2_outputs"
    )
    standAloneProcessDoc1("scatter-valuefrom-wf2") should not equal standAloneProcessDoc2(
        "scatter-valuefrom-wf2"
    )
  }

  "ProcessTranslator" should "render different cwl code for subworkflow but not unchanged tasks" in {

    val standAloneProcessDoc1 = getProcessDocContent("count-lines8-wf")
    val standAloneProcessDoc2 = getProcessDocContent("count-lines8-wf-dup")
    standAloneProcessDoc1("count-lines1-wf.cwl@step_step1@wc-tool.cwl") should equal(
        standAloneProcessDoc2("count-lines1-wf.cwl@step_step1@wc-tool.cwl")
    )
    standAloneProcessDoc1("count-lines1-wf.cwl@step_step2@parseInt-tool.cwl") should equal(
        standAloneProcessDoc2("count-lines1-wf.cwl@step_step2@parseInt-tool.cwl")
    )
    standAloneProcessDoc1("count-lines8-wf.cwl@step_step1@count-lines1-wf.cwl") should equal(
        standAloneProcessDoc2("count-lines8-wf.cwl@step_step1@count-lines1-wf.cwl")
    )

    // step name changed
    standAloneProcessDoc1("count-lines8-wf.cwl_frag_stage-0") should not equal (
        standAloneProcessDoc2("count-lines8-wf.cwl_frag_stage-0")
    )
    standAloneProcessDoc1("count-lines8-wf.cwl_common") should not equal standAloneProcessDoc2(
        "count-lines8-wf.cwl_common"
    )
    standAloneProcessDoc1("count-lines8-wf.cwl_outputs") should not equal standAloneProcessDoc2(
        "count-lines8-wf.cwl_outputs"
    )
    standAloneProcessDoc1("count-lines8-wf.cwl") should not equal standAloneProcessDoc2(
        "count-lines8-wf.cwl"
    )
  }
}
