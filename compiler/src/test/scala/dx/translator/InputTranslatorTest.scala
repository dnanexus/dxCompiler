package dx.translator

import java.nio.file.{Path, Paths}
import dx.Assumptions.isLoggedIn
import dx.Tags.EdgeTest
import dx.api._
import dxCompiler.Main.SuccessfulCompileIR
import dx.core.CliUtils.Failure
import dx.core.Constants
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
    Paths.get(getClass.getResource(s"/${dirname}/${basename}").getPath)
  }

  private val cFlags =
    List("--compileMode", "ir", "--quiet", "--project", DxApi.get.currentProjectId.get)

  // make sure we are logged in

  it should "handle one task and two inputs" in {
    val sourceCode = pathFromBasename("input_file", "add.wdl")
    val inputs = pathFromBasename("input_file", "add_inputs.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "deal with a locked workflow" in {
    val sourceCode = pathFromBasename("input_file", "math.wdl")
    val inputs = pathFromBasename("input_file", "math_inputs.json")
    val args = List(sourceCode.toString,
                    "--inputs",
                    inputs.toString,
                    "--locked"
                    //, "--verbose", "--verboseKey", "GenerateIR"
    ) ++ cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "not compile for several applets without a workflow" in {
    val sourceCode = pathFromBasename("input_file", "several_tasks.wdl")
    val inputs = pathFromBasename("input_file", "several_tasks_inputs.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("cannot generate one input file for 2 tasks")
    }
  }

  it should "one input too many" in {
    val sourceCode = pathFromBasename("input_file", "math.wdl")
    val inputs = pathFromBasename("input_file", "math_inputs2.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString, "--locked") ++ cFlags
    val retval = Main.compile(args.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("Could not map all input fields")
    }
  }

  it should "build defaults into applet underneath workflow" in {
    val sourceCode = pathFromBasename("input_file", "population.wdl")
    val defaults = pathFromBasename("input_file", "population_inputs.json")
    val args = List(sourceCode.toString, "-defaults", defaults.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle inputs specified in the json file, but missing in the workflow" in {
    val sourceCode = pathFromBasename("input_file", "missing_args.wdl")
    val inputs = pathFromBasename("input_file", "missing_args_inputs.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]

    // inputs as defaults
    val args2 = List(sourceCode.toString, "-defaults", inputs.toString) ++ cFlags
    Main.compile(args2.toVector) shouldBe a[SuccessfulCompileIR]

    // Input to an applet.
    // Missing argument in a locked workflow should throw an exception.
    val args3 = List(sourceCode.toString, "--inputs", inputs.toString, "--locked") ++ cFlags
    val retval = Main.compile(args3.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("Could not map all input fields")
    }

    // Missing arguments are legal in an unlocked workflow
    val args4 = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    Main.compile(args4.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "support struct inputs" in {
    val sourceCode = pathFromBasename("struct", "Person.wdl")
    val inputs = pathFromBasename("struct", "Person_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "support array of pairs" taggedAs EdgeTest in {
    val sourceCode = pathFromBasename("input_file", "echo_pairs.wdl")
    val inputs = pathFromBasename("input_file", "echo_pairs_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    //        ++ List("--verbose", "--verboseKey", "GenerateIR")
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle an array of structs" in {
    val sourceCode = pathFromBasename("struct", "array_of_structs.wdl")
    val inputs = pathFromBasename("struct", "array_of_structs_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "override default values in input file" in {
    val sourceCode = pathFromBasename("input_file", "override.wdl")
    val inputs = pathFromBasename("input_file", "override_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "WDL map input" in {
    val sourceCode = pathFromBasename("input_file", "map_argument.wdl")
    val inputs = pathFromBasename("input_file", "map_argument_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

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
    val sourceCode = pathFromBasename("input_file", "map_argument2.wdl")
    val inputs = pathFromBasename("input_file", "map_argument2_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

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
        "gnomad" -> "file-GJ6vK700yzZqxJx94kY0q0gz"
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
    val sourceCode = pathFromBasename("input_file", "no_file_key.wdl")
    val inputs = pathFromBasename("input_file", "no_file_key_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "translate manifest inputs" in {
    val sourceCode = pathFromBasename("manifest", "simple_manifest.wdl")
    val inputs = pathFromBasename("manifest", "simple_manifest_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString, "--useManifests") ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle input files with metadata" in {
    val cwlCode = pathFromBasename("input_file", "input-file-metadata.cwl.json")
    val inputs = pathFromBasename("input_file", "input-file-metadata_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val metadataAttr =
      jsInputs.asJsObject.fields("file1").asJsObject.fields("___").asJsObject.fields("metadata")
    metadataAttr shouldBe JsObject(
        "key1" -> JsNumber(1),
        "key2" -> JsString("value")
    )
  }

  it should "handle parameter name with dot" in {
    val cwlCode = pathFromBasename("input_file", "bwa-mem-tool.cwl.json")
    val inputs = pathFromBasename("input_file", "bwa-mem-tool_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields("reference") shouldBe JsObject(
        "___" -> JsObject(
            "checksum" -> JsString("sha1$hash"),
            "size" -> JsNumber(123),
            "type" -> JsString("File"),
            "uri" -> JsObject(
                "$dnanexus_link" -> JsObject(
                    "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq"),
                    "id" -> JsString("file-GJ6xPPj0yzZxXQFP4kf3FkJq")
                )
            )
        )
    )
    fields("reads") shouldBe JsObject(
        "___" -> JsArray(
            JsObject(
                "type" -> JsString("File"),
                "uri" -> JsObject(
                    "$dnanexus_link" -> JsObject(
                        "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq"),
                        "id" -> JsString("file-GJ6xPPj0yzZvpF8Z4kq03qv3")
                    )
                )
            ),
            JsObject(
                "type" -> JsString("File"),
                "uri" -> JsObject(
                    "$dnanexus_link" -> JsObject(
                        "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq"),
                        "id" -> JsString("file-GJ6xPPj0yzZv8Zvk4kzF5z3x")
                    )
                )
            )
        )
    )
    fields("min_std_max_min") shouldBe JsArray(
        JsNumber(1),
        JsNumber(2),
        JsNumber(3),
        JsNumber(4)
    )
    fields("minimum_seed_length") shouldBe JsNumber(3)
    fields("args__46__py") shouldBe JsObject(
        "___" -> JsObject(
            "type" -> JsString("File"),
            "uri" -> JsObject(
                "$dnanexus_link" -> JsObject(
                    "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq"),
                    "id" -> JsString("file-GJ6xQ600yzZxXQFP4kf3FkJz")
                )
            )
        )
    )
  }

  it should "translate cwl inputs" in {
    val cwlCode = pathFromBasename("input_file", "initialwork-path.cwl.json")
    val inputs = pathFromBasename("input_file", "initialwork-path_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "translate cwl inputs for packed file with two tools" in {
    val cwlCode = pathFromBasename("input_file", "echo-tool-packed.cwl.json")
    val inputs = pathFromBasename("input_file", "echo-tool-packed_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile).asJsObject.fields
    jsInputs.size shouldBe 2
    jsInputs("in") shouldBe JsObject("___" -> JsString("hello test env"))
  }

  it should "translate cwl file inputs" in {
    val cwlCode = pathFromBasename("input_file", "io-any-wf.cwl.json")
    val inputs = pathFromBasename("input_file", "io-any-wf_input3.json")
    val args = List(cwlCode.toString, "--locked", "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile).asJsObject.fields
    jsInputs.size shouldBe 2
    jsInputs.keySet shouldBe Set("bar", "bar___dxfiles")
    val fileObj =
      jsInputs("bar").asJsObject.fields("___").asJsObject.fields("wrapped___").asJsObject.fields
    fileObj("type") shouldBe JsString("File")
    fileObj("uri") match {
      case JsObject(fields) if fields.contains("$dnanexus_link") =>
        fields("$dnanexus_link") shouldBe JsObject(
            "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq"),
            "id" -> JsString("file-GJ6xXF00yzZQpfB04kJGyGj0")
        )
      case other =>
        throw new Exception(s"expected dx link not ${other}")
    }
  }

  it should "translate cwl inputs with default file value" in {
    val cwlCode = pathFromBasename("input_file", "bool-empty-inputbinding.cwl.json")
    val inputs = pathFromBasename("input_file", "bool-empty-inputbinding_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "translate cwl inputs for workflow with no steps" in {
    val cwlCode = pathFromBasename("cwl", "any-type-compat.cwl.json")
    val inputs = pathFromBasename("input_file", "any-type-compat_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    // the input types are all 'Any' so they should be treated as objects
    fields.keySet shouldBe Set("input1",
                               "input1___dxfiles",
                               "input2",
                               "input2___dxfiles",
                               "input3",
                               "input3___dxfiles")
  }

  it should "translate cwl directory input" in {
    val cwlCode = pathFromBasename("input_file", "dir.cwl.json")
    val inputs = pathFromBasename("input_file", "dir_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "translate cwl directory listing I" in {
    val cwlCode = pathFromBasename("input_file", "cat-from-dir.cwl.json")
    val inputs = pathFromBasename("input_file", "cat-from-dir_input1.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields("dir1") shouldBe JsObject(
        Constants.ComplexValueKey -> JsObject(
            "type" -> JsString("Listing"),
            "basename" -> JsString("cwl"),
            "listing" -> JsArray(
                JsObject(
                    "type" -> JsString("File"),
                    "uri" -> JsObject(
                        "$dnanexus_link" -> JsObject(
                            "id" -> JsString("file-GJ6xXF00yzZQpfB04kJGyGj0"),
                            "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")
                        )
                    )
                )
            )
        )
    )
  }

  it should "translate cwl directory listing II" in {
    val cwlCode = pathFromBasename("input_file", "cat-from-dir.cwl.json")
    val inputs = pathFromBasename("input_file", "cat-from-dir_input2.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields("dir1") shouldBe JsObject(
        Constants.ComplexValueKey -> JsObject(
            "type" -> JsString("Listing"),
            "basename" -> JsString("cwl"),
            "listing" -> JsArray(
                JsObject(
                    "type" -> JsString("File"),
                    "uri" -> JsString("literal.txt"),
                    "basename" -> JsString("literal.txt"),
                    "contents" -> JsString("I'm a File literal; howdy!")
                )
            )
        )
    )
  }

  it should "translate cwl string input into enum" in {
    val cwlCode = pathFromBasename("input_file", "enum-string.cwl.json")
    val inputs = pathFromBasename("input_file", "enum-string_input.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields.size shouldBe 1
    fields.keySet should contain("stage-common.example_in")
  }

  it should "translate cwl file with secondary files" in {
    val cwlCode = pathFromBasename("input_file", "dir4.cwl.json")
    val inputs = pathFromBasename("input_file", "dir4_input1.json")
    val args = List(cwlCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields("inf") shouldBe JsObject(
        Constants.ComplexValueKey -> JsObject(
            "type" -> JsString("File"),
            "uri" -> JsObject(
                "$dnanexus_link" -> JsObject(
                    "id" -> JsString("file-GJ6xZ0Q0yzZxXQFP4kf3FkKV"),
                    "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")
                )
            ),
            "secondaryFiles" -> JsArray(
                JsObject(
                    "type" -> JsString("File"),
                    "uri" -> JsObject(
                        "$dnanexus_link" -> JsObject(
                            "id" -> JsString("file-GJ6xg0Q0yzZqxJx94kY0q0p0"),
                            "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")
                        )
                    )
                ),
                JsObject(
                    "type" -> JsString("Folder"),
                    "basename" -> JsString("xtestdir"),
                    "uri" -> JsString("dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:/test_data/cwl/test/")
                )
            )
        )
    )
  }

  it should "translate a WDL file with overrides" in {
    val sourceCode = pathFromBasename("input_file", "runtime_override.wdl")
    val inputs = pathFromBasename("input_file", "runtime_override_inputs.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields.size shouldBe 2
    fields("s") shouldBe JsString("foo")
    fields(Constants.Overrides.encoded) shouldBe JsObject(
        "___" -> JsObject(
            "runtime" -> JsObject(
                "docker" -> JsString("debian:latest"),
                "cpu" -> JsNumber(8)
            )
        )
    )
  }

  it should "translate a CWL file with overrides" in {
    val sourceCode = pathFromBasename("input_file", "env-tool.cwl.json")
    val inputs = pathFromBasename("input_file", "env-tool_input.yaml")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields.size shouldBe 2
    fields("in") shouldBe JsString("hello test env")
    fields(Constants.Overrides.encoded) shouldBe JsObject(
        "___" -> JsObject(
            "requirements" -> JsArray(
                JsObject(
                    "class" -> JsString("EnvVarRequirement"),
                    "envDef" -> JsArray(
                        JsObject(
                            "envName" -> JsString("TEST_ENV"),
                            "envValue" -> JsString("$(inputs.in)")
                        )
                    )
                )
            )
        )
    )
  }
}
