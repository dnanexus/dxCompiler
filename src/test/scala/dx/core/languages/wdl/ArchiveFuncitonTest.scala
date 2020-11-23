package dx.core.languages.wdl

import java.nio.file.Files

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.eval.{DefaultEvalPaths, Eval, WdlValues}
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import wdlTools.syntax.{SourceLocation, WdlVersion}

class ArchiveFuncitonTest extends AnyFlatSpec with Matchers {
  it should "evaluate archive function" in {
    val inputType1 = WdlTypes.T_Array(WdlTypes.T_File, nonEmpty = false)
    val inputType2 = WdlTypes.T_Optional(WdlTypes.T_Boolean)
    val prototype =
      WdlTypes.T_Function2("archive_file_array", inputType1, inputType2, WdlTypes.T_File)
    ArchiveFunction
      .getPrototype("archive_file_array", Vector(inputType1, inputType2)) shouldBe prototype
    val tmpdir = Files.createTempDirectory("test")
    val file1 = tmpdir.

    val paths = DefaultEvalPaths.createFromTemp()
    val evaluator = Eval(paths, Some(WdlVersion.V1), Vector(ArchiveFunction))
    val expr = TAT.ExprApply(
        "archive_file_array",
        prototype,
        Vector(
            TAT.ExprArray(Vector(TAT.ValueFile(), TAT.ValueFile()),
                          WdlTypes.T_Array(WdlTypes.T_File, nonEmpty = false),
                          SourceLocation.empty)
        ),
        WdlTypes.T_File,
        SourceLocation.empty
    )
    val value = evaluator.applyExpr(expr)
    val path = value match {
      case WdlValues.V_File(path) => path
      case _ => throw new AssertionError("return value is not a V_File")
    }
  }
}
