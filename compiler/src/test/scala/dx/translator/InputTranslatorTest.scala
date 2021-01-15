package dx.translator

import java.nio.file.{Path, Paths}
import dx.Assumptions.isLoggedIn
import dx.Tags.EdgeTest
import dx.api._
import dxCompiler.Main.SuccessIR
import dx.core.CliUtils.Failure
import dx.util.{FileUtils, JsUtils}
import dxCompiler.Main
import org.scalatest.Inside._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

// These tests involve compilation -without- access to the platform.
//
class InputTranslatorTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)

  private def pathFromBasename(dirname: String, basename: String): Path = {
    val p = getClass.getResource(s"/${dirname}/${basename}").getPath
    Paths.get(p)
  }

  private val cFlags =
    List("--compileMode", "ir", "-quiet", "--project", DxApi.get.currentProject.id)

  // make sure we are logged in

  it should "handle one task and two inputs" in {
    val wdlCode = pathFromBasename("input_file", "add.wdl")
    val inputs = pathFromBasename("input_file", "add_inputs.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    Main.compile(args.toVector) shouldBe a[SuccessIR]
  }

  it should "deal with a locked workflow" in {
    val wdlCode = pathFromBasename("input_file", "math.wdl")
    val inputs = pathFromBasename("input_file", "math_inputs.json")
    val args = List(wdlCode.toString,
                    "-inputs",
                    inputs.toString,
                    "--locked"
                    //, "--verbose", "--verboseKey", "GenerateIR"
    ) ++ cFlags
    Main.compile(args.toVector) shouldBe a[SuccessIR]
  }

  it should "not compile for several applets without a workflow" in {
    val wdlCode = pathFromBasename("input_file", "several_tasks.wdl")
    val inputs = pathFromBasename("input_file", "several_tasks_inputs.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("Cannot generate one input file for 2 tasks")
    }
  }

  it should "one input too many" in {
    val wdlCode = pathFromBasename("input_file", "math.wdl")
    val inputs = pathFromBasename("input_file", "math_inputs2.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString, "--locked") ++ cFlags
    val retval = Main.compile(args.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("Could not map all input fields")
    }
  }

  it should "build defaults into applet underneath workflow" in {
    val wdlCode = pathFromBasename("input_file", "population.wdl")
    val defaults = pathFromBasename("input_file", "population_inputs.json")
    val args = List(wdlCode.toString, "-defaults", defaults.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]
  }

  it should "handle inputs specified in the json file, but missing in the workflow" in {
    val wdlCode = pathFromBasename("input_file", "missing_args.wdl")
    val inputs = pathFromBasename("input_file", "missing_args_inputs.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    Main.compile(args.toVector) shouldBe a[SuccessIR]

    // inputs as defaults
    val args2 = List(wdlCode.toString, "-defaults", inputs.toString) ++ cFlags
    Main.compile(args2.toVector) shouldBe a[SuccessIR]

    // Input to an applet.
    // Missing argument in a locked workflow should throw an exception.
    val args3 = List(wdlCode.toString, "-inputs", inputs.toString, "--locked") ++ cFlags
    val retval = Main.compile(args3.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("Could not map all input fields")
    }

    // Missing arguments are legal in an unlocked workflow
    val args4 = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    Main.compile(args4.toVector) shouldBe a[SuccessIR]
  }

  it should "support struct inputs" in {
    val wdlCode = pathFromBasename("struct", "Person.wdl")
    val inputs = pathFromBasename("struct", "Person_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]
  }

  it should "support array of pairs" taggedAs EdgeTest in {
    val wdlCode = pathFromBasename("input_file", "echo_pairs.wdl")
    val inputs = pathFromBasename("input_file", "echo_pairs_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    //        ++ List("--verbose", "--verboseKey", "GenerateIR")
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]
  }

  it should "handle an array of structs" in {
    val wdlCode = pathFromBasename("struct", "array_of_structs.wdl")
    val inputs = pathFromBasename("struct", "array_of_structs_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]
  }

  it should "override default values in input file" in {
    val wdlCode = pathFromBasename("input_file", "override.wdl")
    val inputs = pathFromBasename("input_file", "override_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]
  }

  it should "WDL map input" in {
    val wdlCode = pathFromBasename("input_file", "map_argument.wdl")
    val inputs = pathFromBasename("input_file", "map_argument_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields.size shouldBe 2 // there's an extra input for the ___files array
    fields.keySet should contain("stage-common.names_and_phones")
    val mapObj = fields("stage-common.names_and_phones") match {
      case JsObject(fields2) =>
        fields2.get("___") match {
          case Some(JsObject(fields3)) => fields3
          case _                       => throw new AssertionError("expected ___ object")
        }
      case _ => throw new AssertionError("expected object")
    }
    mapObj.keySet shouldBe Set("keys", "values")
    val JsArray(keys) = mapObj("keys")
    val JsArray(values) = mapObj("values")
    keys.size shouldEqual values.size
    val map = keys.zip(values).toMap
    map should contain theSameElementsAs Map(
        JsString("Gina") -> JsNumber(334320),
        JsString("Chad") -> JsNumber(578333),
        JsString("Valery") -> JsNumber(88339)
    )
  }

  it should "WDL map input with file values" in {
    val wdlCode = pathFromBasename("input_file", "map_argument2.wdl")
    val inputs = pathFromBasename("input_file", "map_argument2_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields.size shouldBe 2 // there's an extra input for the ___files array
    fields.keySet should contain("stage-common.files")
    val mapObj = fields("stage-common.files") match {
      case JsObject(fields2) =>
        fields2.get("___") match {
          case Some(JsObject(fields3)) => fields3
          case _                       => throw new AssertionError("expected ___ object")
        }
      case _ => throw new AssertionError("expected object")
    }
    mapObj.keySet shouldBe Set("keys", "values")
    val JsArray(keys) = mapObj("keys")
    val JsArray(values) = mapObj("values")
    keys.size shouldEqual values.size
    val map = keys.zip(values).toMap
    val fileIds = Map(
        "splicing" -> "file-FV5fqXj0ffPB9bKP986j5kVQ",
        "transcripts" -> "file-Fg5PgBQ0ffP7B8bg3xqB115G",
        "cosmic" -> "file-Fg5PgBj0ffPP0Jjv3zfv0yxq",
        "dbsnp" -> "file-Fy9VJ1j0yzZgVgFqJPf6KK17",
        "evs" -> "file-FGqFJ8Q0ffPGVz3zGy4FK02P",
        "cmh_maf" -> "file-FGzzpkQ0ffPJX74548Vp6670",
        "exac" -> "file-FGqFY200ffP3qQYgK163Z4gf",
        "gnomad" -> "file-Fzy8x5j0yzZVQv9KB41FGz3V"
    )
    map.foreach {
      case (JsString(key), JsObject(value)) =>
        val fileId = value("$dnanexus_link") match {
          case JsString(fileId) => fileId
          case JsObject(fields) =>
            val JsString(fileId) = fields("id")
            fileId
          case _ => throw new AssertionError("invalid file value")
        }
        fileId shouldBe fileIds(key)
      case _ => throw new AssertionError("expected string -> object")
    }
  }

  it should "allow file as WDL map key" in {
    val wdlCode = pathFromBasename("input_file", "no_file_key.wdl")
    val inputs = pathFromBasename("input_file", "no_file_key_input.json")
    val args = List(wdlCode.toString, "-inputs", inputs.toString, "-verbose") ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessIR]
  }
}
