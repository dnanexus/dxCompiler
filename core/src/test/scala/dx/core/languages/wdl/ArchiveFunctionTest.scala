package dx.core.languages.wdl

import java.nio.file.{Files, Paths}

import dx.core.ir.{LocalizedArchive, PackedArchive, Type, Value}
import dx.util.FileUtils
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.eval.{DefaultEvalPaths, Eval, WdlValues}
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import wdlTools.syntax.{SourceLocation, WdlVersion}

class ArchiveFunctionTest extends AnyFlatSpec with Matchers {
  it should "evaluate archive function" in {
    val inputType1 = WdlTypes.T_Array(WdlTypes.T_File, nonEmpty = false)
    val inputType2 = WdlTypes.T_Optional(WdlTypes.T_Boolean)
    val prototype =
      WdlTypes.T_Function2("archive_file_array", inputType1, inputType2, WdlTypes.T_File)
    ArchiveFunction
      .getPrototype("archive_file_array", Vector(inputType1, inputType2)) shouldBe Some(prototype)

    val paths = DefaultEvalPaths.createFromTemp()
    val file1 = paths.getRootDir().resolve("file1.txt")
    FileUtils.writeFileContent(file1, "file1")
    val file2 = paths.getRootDir().resolve("file2.txt")
    FileUtils.writeFileContent(file2, "file2")

    val evaluator = Eval(paths, Some(WdlVersion.V1), Vector(ArchiveFunction))
    val expr = TAT.ExprApply(
        "archive_file_array",
        prototype,
        Vector(
            TAT.ExprArray(
                Vector(TAT.ValueFile(file1.toString, WdlTypes.T_File, SourceLocation.empty),
                       TAT.ValueFile(file2.toString, WdlTypes.T_File, SourceLocation.empty)),
                WdlTypes.T_Array(WdlTypes.T_File, nonEmpty = false),
                SourceLocation.empty
            ),
            TAT.ValueBoolean(value = true, WdlTypes.T_Boolean, SourceLocation.empty)
        ),
        WdlTypes.T_File,
        SourceLocation.empty
    )
    val value = evaluator.applyExpr(expr)
    val path = value match {
      case WdlValues.V_File(path) => path
      case _                      => throw new AssertionError("return value is not a V_File")
    }

    Files.exists(file1) shouldBe false
    Files.exists(file2) shouldBe false

    val packed = PackedArchive(Paths.get(path))()
    val (localized, _) = packed.localize()
    try {
      localized.irValue match {
        case Value.VArray(Vector(file1: Value.VFile, file2: Value.VFile)) =>
          FileUtils.readFileContent(Paths.get(file1.uri)) shouldBe "file1"
          FileUtils.readFileContent(Paths.get(file2.uri)) shouldBe "file2"
        case _ =>
          throw new Exception(s"unexpected IR value ${localized.irValue}")
      }
    } finally {
      localized.close()
    }
  }

  it should "evaluate unarchive function" in {
    val paths = DefaultEvalPaths.createFromTemp()
    val file1 = paths.getRootDir(ensureExists = true).resolve("file1.txt")
    FileUtils.writeFileContent(file1, "file1")
    val file2 = paths.getRootDir(ensureExists = true).resolve("file2.txt")
    FileUtils.writeFileContent(file2, "file2")

    val irType = Type.TArray(Type.TFile)
    val irValue = Value.VArray(Vector(Value.VFile(file1.toString), Value.VFile(file2.toString)))
    val archive =
      LocalizedArchive(irType, irValue)(parentDir = Some(paths.getRootDir(ensureExists = true)))
    val packedArchive = archive.pack()

    val prototype =
      WdlTypes.T_Function1("unarchive_file_array", WdlTypes.T_File, WdlTypes.T_Any)
    UnarchiveFunction
      .getPrototype("unarchive_file_array", Vector(WdlTypes.T_File)) shouldBe Some(prototype)

    val evaluator = Eval(paths, Some(WdlVersion.V1), Vector(UnarchiveFunction))
    val expr = TAT.ExprApply(
        "unarchive_file_array",
        prototype,
        Vector(
            TAT.ValueFile(packedArchive.path.toString, WdlTypes.T_File, SourceLocation.empty)
        ),
        WdlTypes.T_Any,
        SourceLocation.empty
    )
    evaluator.applyExpr(expr) match {
      case WdlValues.V_Array(Vector(WdlValues.V_File(path1), WdlValues.V_File(path2))) =>
        FileUtils.readFileContent(Paths.get(path1)) shouldBe "file1"
        FileUtils.readFileContent(Paths.get(path2)) shouldBe "file2"
      case other =>
        throw new AssertionError(s"unexpected expression value ${other}")
    }
  }
}
