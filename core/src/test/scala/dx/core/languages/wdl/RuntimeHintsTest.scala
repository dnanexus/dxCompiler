package dx.core.languages.wdl

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.eval.{DefaultEvalPaths, Eval, EvalException}
import wdlTools.eval.WdlValues._
import wdlTools.syntax.{Quoting, SourceLocation, WdlVersion}
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import dx.util.{FileSourceResolver, Logger}

import scala.collection.immutable.SeqMap

class RuntimeHintsTest extends AnyFlatSpec with Matchers {
  private val evaluator: Eval =
    Eval(DefaultEvalPaths.empty,
         Some(WdlVersion.V1),
         Vector.empty,
         FileSourceResolver.get,
         Logger.get)

  private def stringExpr(s: String): TAT.Expr =
    TAT.ValueString(s, WdlTypes.T_String, Quoting.Double)(SourceLocation.empty)

  private def runtimeWith(entries: (String, String)*): Runtime = {
    val rt = SeqMap(entries.map { case (k, v) => k -> stringExpr(v) }: _*)
    Runtime(
        WdlVersion.V1,
        Some(TAT.RuntimeSection(rt)(SourceLocation.empty)),
        None,
        evaluator
    )
  }

  private def hintsRuntimeV2(dxFields: (String, String)*): Runtime = {
    val dnanexusInner = SeqMap(dxFields.map {
      case (k, v) => k -> TAT.MetaValueString(v, Quoting.Double)(SourceLocation.empty)
    }: _*)
    val hints = SeqMap(
        Runtime.DxHintsKey ->
          TAT.MetaValueObject(dnanexusInner)(SourceLocation.empty)
    )
    Runtime(
        WdlVersion.V2,
        None,
        Some(TAT.MetaSection(hints)(SourceLocation.empty)),
        evaluator
    )
  }

  it should "return None when dx_shm_size is not set" in {
    runtimeWith().shmSize shouldBe None
    runtimeWith().ipcMode shouldBe None
  }

  it should "extract dx_shm_size from runtime block (WDL 1.x)" in {
    runtimeWith(Runtime.DxShmSizeKey -> "8g").shmSize shouldBe Some("8g")
  }

  it should "extract dx_ipc_mode from runtime block (WDL 1.x)" in {
    runtimeWith(Runtime.DxIpcModeKey -> "host").ipcMode shouldBe Some("host")
  }

  it should "extract shm_size and ipc_mode from hints.dnanexus block (WDL 2.0)" in {
    val rt = hintsRuntimeV2("shm_size" -> "4g", "ipc_mode" -> "shareable")
    rt.shmSize shouldBe Some("4g")
    rt.ipcMode shouldBe Some("shareable")
  }

  it should "reject malformed dx_shm_size values" in {
    val rt = runtimeWith(Runtime.DxShmSizeKey -> "garbage")
    val ex = intercept[EvalException](rt.shmSize)
    ex.getMessage should include(Runtime.DxShmSizeKey)
  }

  it should "reject malformed dx_ipc_mode values" in {
    val rt = runtimeWith(Runtime.DxIpcModeKey -> "weird-value")
    val ex = intercept[EvalException](rt.ipcMode)
    ex.getMessage should include(Runtime.DxIpcModeKey)
  }

  it should "accept all valid shm_size suffixes" in {
    Seq("64", "64b", "64k", "64m", "8g", "1024M", "2G").foreach { v =>
      noException should be thrownBy Runtime.validateShmSize(v)
    }
  }

  it should "reject zero and leading-zero shm_size values" in {
    Seq("0", "0g", "00", "0064m").foreach { v =>
      an[EvalException] should be thrownBy Runtime.validateShmSize(v)
    }
  }

  it should "accept all valid ipc_mode values" in {
    Seq("none", "private", "shareable", "host", "container:my-container").foreach { v =>
      noException should be thrownBy Runtime.validateIpcMode(v)
    }
  }

  it should "reject ipc_mode values that could enable shell injection" in {
    // The container:.+ form is the security boundary: anything beyond Docker's container-name
    // grammar (https://docs.docker.com/reference/cli/docker/container/run/#name) could allow
    // a malicious WDL author to inject extra docker flags via the shell-rendered run command.
    Seq(
        "container:foo --privileged",
        "container:foo;rm -rf /",
        "container:foo$(whoami)",
        "container:foo`id`",
        "container:",
        "container:.bad",
        "weird-value",
        "host;rm -rf /"
    ).foreach { v =>
      an[EvalException] should be thrownBy Runtime.validateIpcMode(v)
    }
  }

  it should "exercise the WDL 1.0 customer scenario from APPS-3954" in {
    // Customer's WDL has `runtime { dx_shm_size: "8g" }` directly, no hints block.
    // This trace pins down: V1 runtime block -> getDxHint(ShmSize) -> Some("8g").
    val rt = runtimeWith(Runtime.DxShmSizeKey -> "8g")
    rt.shmSize shouldBe Some("8g")
    rt.getDxHint(Runtime.ShmSize) shouldBe Some(V_String("8g"))
  }
}
