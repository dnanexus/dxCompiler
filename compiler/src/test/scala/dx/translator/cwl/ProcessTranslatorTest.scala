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
    val deconstructedCallables1 = getDeconstructedProcesses("count-lines1-wf")
    val deconstructedCallables2 = getDeconstructedProcesses("count-lines1-wf-dup")
    deconstructedCallables1("parseInt-tool.cwl") should equal(
        deconstructedCallables2("parseInt-tool.cwl")
    )

    deconstructedCallables1("wc-tool.cwl") should not equal deconstructedCallables2("wc-tool.cwl")
    deconstructedCallables1("count-lines1-wf") should not equal deconstructedCallables2(
        "count-lines1-wf"
    )

  }

  "ProcessTranslator" should "render different cwl code for changed common/outputs" in {
    // TODO: not cover the case when multiple process/step has the same simpilied ids
    // need to fix getDocContent function
    val deconstructedCallables1 = getDeconstructedProcesses("scatter-valuefrom-wf2")
    val deconstructedCallables2 = getDeconstructedProcesses("scatter-valuefrom-wf2-dup")
    deconstructedCallables1("step1command") should equal(deconstructedCallables2("step1command"))
    deconstructedCallables1("scatter-valuefrom-wf2_frag_stage-0") should equal(
        deconstructedCallables2("scatter-valuefrom-wf2_frag_stage-0")
    )

    deconstructedCallables1("scatter-valuefrom-wf2_common") should not equal deconstructedCallables2(
        "scatter-valuefrom-wf2_common"
    )
    deconstructedCallables1("scatter-valuefrom-wf2_outputs") should not equal deconstructedCallables2(
        "scatter-valuefrom-wf2_outputs"
    )
    deconstructedCallables1("scatter-valuefrom-wf2") should not equal deconstructedCallables2(
        "scatter-valuefrom-wf2"
    )
  }

  // it should "render same wdl code for every unchanged block/app/frag and different if there are changes" in {
  //   val (docV1, typeAliasesV1, versionSupportV1) =
  //     VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v1.wdl"))
  //   val (docV2, typeAliasesV2, versionSupportV2) =
  //     VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v2.wdl"))
  //   val deconstructedCallablesV1 = getDeconstructedCallables(
  //       doc = docV1,
  //       typeAliases = typeAliasesV1,
  //       versionSupport = versionSupportV1
  //   )
  //   val deconstructedCallablesV2 = getDeconstructedCallables(
  //       doc = docV2,
  //       typeAliases = typeAliasesV2,
  //       versionSupport = versionSupportV2
  //   )

  //   deconstructedCallablesV1("reuse_print") should equal(deconstructedCallablesV2("reuse_print"))
  //   deconstructedCallablesV1("reuse_block_2") should equal(
  //       deconstructedCallablesV2("reuse_block_2")
  //   )
  //   deconstructedCallablesV1("reuse_frag_stage-12") should equal(
  //       deconstructedCallablesV2(
  //           "reuse_frag_stage-12"
  //       )
  //   )
  //   deconstructedCallablesV1("reuse_frag_stage-0") should not equal (deconstructedCallablesV2(
  //       "reuse_frag_stage-0"
  //   ))
  //   deconstructedCallablesV1("reuse") should not equal (deconstructedCallablesV2(
  //       "reuse"
  //   ))
  // }
}
