package dx.translator.wdl

import dx.core.ir.{Application, Callable, InstanceTypeSelection, Workflow}
import dx.core.languages.wdl.{VersionSupport, WdlBundle}
import dx.translator.DefaultReorgSettings
import dx.util.Bindings
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}

import java.nio.file.{Path, Paths}

class CallableTranslatorTest extends AnyFlatSpec with Matchers {

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  private def getDeconstructedCallables(doc: TAT.Document,
                                        typeAliases: Bindings[String, WdlTypes.T_Struct],
                                        versionSupport: VersionSupport): Map[String, String] = {
    val wdlBundle: WdlBundle = WdlBundle.create(doc = doc)
    val callableTranslator: CallableTranslator = CallableTranslator(
        wdlBundle = wdlBundle,
        typeAliases = typeAliases.toMap,
        locked = false,
        defaultRuntimeAttrs = Map.empty,
        reorgAttrs = DefaultReorgSettings(false),
        perWorkflowAttrs = Map.empty,
        defaultScatterChunkSize = 1,
        useManifests = false,
        instanceTypeSelection = InstanceTypeSelection.Static,
        versionSupport = versionSupport
    )
    val depOrder: Vector[TAT.Callable] = wdlBundle.sortByDependencies()
    val (_, sortedCallables) =
      depOrder.foldLeft((Map.empty[String, Callable], Vector.empty[Callable])) {
        case ((allCallables, sortedCallables), callable) =>
          val translatedCallables = callableTranslator.translateCallable(callable, allCallables)
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

  "CallableTranslator" should "render different wdl code for every block/app/frag" in {
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v1.wdl"))
    val deconstructedCallables = getDeconstructedCallables(
        doc = doc,
        typeAliases = typeAliases,
        versionSupport = versionSupport
    )
    deconstructedCallables("reuse") should equal(deconstructedCallables("reuse"))
    deconstructedCallables("reuse_print") should not equal deconstructedCallables("reuse_multiply")
    deconstructedCallables("reuse_block_2") should not equal deconstructedCallables("reuse_block_4")
    deconstructedCallables("reuse_frag_stage-12") should not equal deconstructedCallables(
        "reuse_frag_stage-6"
    )
    deconstructedCallables("reuse_frag_stage-0") should not equal deconstructedCallables(
        "reuse_block_4"
    )
  }

  it should "render same wdl code for every unchanged block/app/frag and different if there are changes" in {
    val (docV1, typeAliasesV1, versionSupportV1) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v1.wdl"))
    val (docV2, typeAliasesV2, versionSupportV2) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v2.wdl"))
    val deconstructedCallablesV1 = getDeconstructedCallables(
        doc = docV1,
        typeAliases = typeAliasesV1,
        versionSupport = versionSupportV1
    )
    val deconstructedCallablesV2 = getDeconstructedCallables(
        doc = docV2,
        typeAliases = typeAliasesV2,
        versionSupport = versionSupportV2
    )

    deconstructedCallablesV1("reuse_print") should equal(deconstructedCallablesV2("reuse_print"))
    deconstructedCallablesV1("reuse_block_2") should equal(
        deconstructedCallablesV2("reuse_block_2")
    )
    deconstructedCallablesV1("reuse_frag_stage-12") should equal(
        deconstructedCallablesV2(
            "reuse_frag_stage-12"
        )
    )
    deconstructedCallablesV1("reuse_frag_stage-0") should not equal (deconstructedCallablesV2(
        "reuse_frag_stage-0"
    ))
    deconstructedCallablesV1("reuse") should not equal (deconstructedCallablesV2(
        "reuse"
    ))
  }
}
