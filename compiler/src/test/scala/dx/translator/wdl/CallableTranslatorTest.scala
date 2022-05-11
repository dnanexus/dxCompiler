package dx.translator.wdl

import dx.api.InstanceTypeRequest
import dx.core.ir.RunSpec.StaticInstanceType
import dx.core.ir.{Application, Callable, InstanceTypeSelection, Parameter, Type, Workflow}
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

  private def getSortedCallables(doc: TAT.Document,
                                 typeAliases: Bindings[String, WdlTypes.T_Struct],
                                 versionSupport: VersionSupport): Vector[Callable] = {
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
    sortedCallables
  }

  private def deconstructCallables(sortedCallables: Vector[Callable]): Map[String, String] = {
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
    val deconstructedCallables = deconstructCallables(
        getSortedCallables(
            doc = doc,
            typeAliases = typeAliases,
            versionSupport = versionSupport
        )
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

  // APPS-1175 - frag
  it should "translate with a workflow wrapped in a frag applet and propagate the outputs of its stages to the " +
    "outputs of the frag" in {
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceFile(
          pathFromBasename("subworkflows", "apps_1175_nested_wf_frag.wdl")
      )
    val sortedCallables = getSortedCallables(
        doc = doc,
        typeAliases = typeAliases,
        versionSupport = versionSupport
    )
    sortedCallables.count(_.name.contains("frag")) shouldBe 1
    val fragWrapperOuts = sortedCallables.filter(_.name.contains("frag")).head match {
      case app: Application => app.outputs
      case wf: Workflow     => wf.outputs
      case _                => throw new Exception("Unexpected callable type")
    }
    fragWrapperOuts.size shouldBe 2
    fragWrapperOuts shouldBe a[Vector[_]]
    val outTypes: Map[String, Type] = (fragWrapperOuts map {
      case (param: Parameter, _) => param.name.toString -> param.dxType
      case param: Parameter      => param.name.toString -> param.dxType
    }).toMap
    outTypes("nested_inner.nested_inner_wf_out") shouldBe Type.TFile
    outTypes("nested_inner.___test_inner1.test_out") shouldBe Type.TFile
  }

  // APPS-1175 - frag with 2 separate calls instead of a nested workflow
  it should "translate a frag wrapper with 2 calls and maintain outputs from both calls" in {
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceFile(
          pathFromBasename("subworkflows", "apps_1175_frag_2_calls.wdl")
      )
    val sortedCallables = getSortedCallables(
        doc = doc,
        typeAliases = typeAliases,
        versionSupport = versionSupport
    )
    sortedCallables.count(_.name.contains("frag")) shouldBe 1
    val fragWrapperOuts = sortedCallables.filter(_.name.contains("frag")).head match {
      case app: Application => app.outputs
      case wf: Workflow     => wf.outputs
      case _                => throw new Exception("Unexpected callable type")
    }
    fragWrapperOuts.size shouldBe 3
    val outTypes: Map[String, Type] = (fragWrapperOuts map {
      case (param: Parameter, _) => param.name.toString -> param.dxType
      case param: Parameter      => param.name.toString -> param.dxType
    }).toMap
    outTypes("nested_inner.nested_inner_wf_out") shouldBe Type.TFile
    outTypes("nested_inner.___test_inner1.test_out") shouldBe Type.TFile
  }

  // APPS-1175 - conditional
  it should "translate with a workflow wrapped in a conditional frag applet and propagate the outputs of its stages " +
    "to the outputs of the frag" in {
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceFile(
          pathFromBasename("subworkflows", "apps_1175_nested_wf_conditional.wdl")
      )
    val sortedCallables = getSortedCallables(
        doc = doc,
        typeAliases = typeAliases,
        versionSupport = versionSupport
    )
    sortedCallables.count(_.name.contains("frag")) shouldBe 1
    val fragWrapperOuts = sortedCallables.filter(_.name.contains("frag")).head match {
      case app: Application => app.outputs
      case wf: Workflow     => wf.outputs
      case _                => throw new Exception("Unexpected callable type")
    }
    fragWrapperOuts.size shouldBe 2
    fragWrapperOuts shouldBe a[Vector[_]]
    val outTypes: Map[String, Type] = (fragWrapperOuts map {
      case (param: Parameter, _) => param.name.toString -> param.dxType
      case param: Parameter      => param.name.toString -> param.dxType
    }).toMap
    outTypes("nested_inner.nested_inner_wf_out") shouldBe Type.TFile
    outTypes("nested_inner.___test_inner1.test_out") shouldBe Type.TFile
  }

  it should "render same wdl code for every unchanged block/app/frag and different if there are changes" in {
    val (docV1, typeAliasesV1, versionSupportV1) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v1.wdl"))
    val (docV2, typeAliasesV2, versionSupportV2) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_994_v2.wdl"))
    val deconstructedCallablesV1 = deconstructCallables(
        getSortedCallables(
            doc = docV1,
            typeAliases = typeAliasesV1,
            versionSupport = versionSupportV1
        )
    )
    val deconstructedCallablesV2 = deconstructCallables(
        getSortedCallables(
            doc = docV2,
            typeAliases = typeAliasesV2,
            versionSupport = versionSupportV2
        )
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

  it should "render frag wrapper to store the requested instance type of the wrapped child process" in {
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_1128_wrap_native_exec.wdl"))
    val sortedCallables = getSortedCallables(
        doc = doc,
        typeAliases = typeAliases,
        versionSupport = versionSupport
    )
    val wrapperFrag = sortedCallables.filter(_.name.contains("frag_stage")).head
    wrapperFrag match {
      case app: Application =>
        app.instanceType shouldBe StaticInstanceType(
            InstanceTypeRequest(Some("mem1_ssd1_v2_x8"),
                                None,
                                None,
                                None,
                                None,
                                None,
                                None,
                                None,
                                None,
                                None)
        )
      case _ => throw new Exception("Unexpected Value")
    }
  }

}
