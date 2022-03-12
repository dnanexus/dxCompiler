package dx.translator.wdl

import dx.core.ir.{Application, Callable, InstanceTypeSelection, Workflow}
import dx.core.languages.wdl.{VersionSupport, WdlBundle}
import dx.translator.DefaultReorgSettings
import dx.util.CodecUtils
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.types.{TypedAbstractSyntax => TAT}

import java.nio.file.{Path, Paths}

class CallableTranslatorTest extends AnyFlatSpec with Matchers {

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }
  "CallableTranslator" should "render different wdl code for every block/app/frag" in {
    val (docV1, typeAliasesV1, versionSupportV1) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v1.wdl"))
    val wdlBundleV1: WdlBundle = WdlBundle.create(doc = docV1)
    val callableTranslatorV1: CallableTranslator = CallableTranslator(
        wdlBundle = wdlBundleV1,
        typeAliases = typeAliasesV1.toMap,
        locked = false,
        defaultRuntimeAttrs = Map.empty,
        reorgAttrs = DefaultReorgSettings(false),
        perWorkflowAttrs = Map.empty,
        defaultScatterChunkSize = 1,
        useManifests = false,
        instanceTypeSelection = InstanceTypeSelection.Static,
        versionSupport = versionSupportV1
    )
    val depOrder: Vector[TAT.Callable] = wdlBundleV1.sortByDependencies()
    val (_, sortedCallables) =
      depOrder.foldLeft((Map.empty[String, Callable], Vector.empty[Callable])) {
        case ((allCallables, sortedCallables), callable) =>
          val translatedCallables = callableTranslatorV1.translateCallable(callable, allCallables)
          (
              allCallables ++ translatedCallables.map(c => c.name -> c).toMap,
              sortedCallables ++ translatedCallables
          )
      }
    val deconstructedCallables: Map[String, String] = sortedCallables.map {
      case Application(name, _, _, _, _, _, document, _, _, _, _, _) =>
        name -> CodecUtils.md5Checksum(document.toString)
      case Workflow(name, _, _, _, document, _, _, _, _, _, _) =>
        name -> CodecUtils.md5Checksum(document.toString)
    }.toMap
    deconstructedCallables("reuse_print") should not equal deconstructedCallables("reuse_multiply")
    deconstructedCallables("reuse_block_2") should not equal deconstructedCallables("reuse_block_4")
    deconstructedCallables("reuse_frag_stage-12") should not equal deconstructedCallables(
        "reuse_frag_stage-6"
    )
    deconstructedCallables(" reuse_frag_stage-0") should not equal deconstructedCallables(
        "reuse_block_4"
    )
  }
}
