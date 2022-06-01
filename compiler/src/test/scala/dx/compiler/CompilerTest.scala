package dx.compiler

import java.io.File
import java.nio.file.{Path, Paths}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import java.util.UUID.randomUUID
import dx.Assumptions.{isLoggedIn, toolkitCallable}
import dx.Tags.NativeTest
import dx.api._
import dx.core.Constants
import dx.core.ir.Callable
import dx.core.CliUtils.Termination
import dx.util.{FileUtils, JsUtils, Logger, SysUtils}
import dxCompiler.Main
import dxCompiler.Main.{SuccessfulCompileIR, SuccessfulCompileNativeNoTree}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

// This test module requires being logged in to the platform. It compiles WDL scripts without the runtime library.
// This tests the compiler Native mode, however, it creates DNAnexus applets and workflows that are not runnable.
class CompilerTest extends AnyFlatSpec with Matchers with BeforeAndAfterAll {
  assume(isLoggedIn)
  assume(toolkitCallable)
  private val logger = Logger.Quiet
  private val dxApi = DxApi()(logger)

  private def pathFromBasename(dir: String, basename: String): Path = {
    val p = getClass.getResource(s"/${dir}/${basename}").getPath
    Paths.get(p)
  }

  val testProject = "dxCompiler_playground"

  private val dxTestProject =
    try {
      dxApi.resolveProject(testProject)
    } catch {
      case _: Exception =>
        throw new Exception(
            s"""|Could not find project ${testProject}, you probably need to be logged into
                |the platform""".stripMargin
        )
    }

  private val username = dxApi.whoami()
  private val unitTestsPath = {
    val dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm")
    val test_time = dateFormatter.format(LocalDateTime.now)
    s"/unit_tests/${username}/CompilerTest/${test_time}_${randomUUID().toString.substring(24)}"
  }
  private val reorgAppletFolder = s"${unitTestsPath}/reorg_applets"
  private val unitTestsReusePath = s"${unitTestsPath}/reuse"
  private var reorgAppletId: Option[String] = None

  private val cFlagsBase: List[String] = List(
      "-project",
      dxTestProject.id,
      "-quiet",
      "-force"
  )
  private val cFlags: List[String] = cFlagsBase ++ List("-compileMode",
                                                        "NativeWithoutRuntimeAsset",
                                                        "-folder",
                                                        unitTestsPath,
                                                        "-locked")

  private val cFlagsReorgIR: List[String] = cFlagsBase ++
    List("-compileMode", "IR", "-folder", "/reorg_tests")
  private val cFlagsReorgCompile: List[String] = cFlagsBase ++
    List("-compileMode", "NativeWithoutRuntimeAsset", "-folder", "/reorg_tests")
  private lazy val cFlagsReuse: List[String] = cFlagsBase ++
    List("-compileMode",
         "NativeWithoutRuntimeAsset",
         "-folder",
         unitTestsReusePath,
         "-projectWideReuse")

  override def beforeAll(): Unit = {
    // build the directory with the native applets
    dxTestProject.newFolder(reorgAppletFolder, parents = true)
    // building necessary applets before starting the tests
    val topDir = Paths.get(System.getProperty("user.dir"))
    val (_, appletIdJs, _) = SysUtils.execCommand(
        s"dx build $topDir/test/applets/functional_reorg_test --destination ${testProject}:${reorgAppletFolder} --brief"
    )
    val JsString(appletId) = JsUtils.jsFromString(appletIdJs).asJsObject.fields("id")
    reorgAppletId = Some(appletId)
  }

  override def afterAll(): Unit = {
    dxTestProject.removeFolder(unitTestsPath, recurse = true)
  }

  private def compileGetStages(fixturePath: String,
                               compileArgs: List[String]): Vector[DxWorkflowStageDesc] = {
    val args = fixturePath :: compileArgs
    val retval = Main.compile(args.toVector)
    val wfId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(wfId)) => wfId
      case other                                          => throw new Exception(s"expected single workflow not ${other}")
    }
    val stages = dxApi
      .workflow(wfId)
      .describe()
      .stages
      .get
    stages
  }

  private object WithExtras {
    def apply(extrasContent: String)(f: String => Termination): Termination = {
      val tmpExtras = File.createTempFile("reorg-", ".json")
      FileUtils.writeFileContent(tmpExtras.toPath, extrasContent)
      try {
        f(tmpExtras.toString)
      } finally {
        tmpExtras.delete()
      }
    }
  }

  "Compiler" should "compile a workflow with a native app wrapped in frag and override instance based using RAM" in {
    val path = pathFromBasename("bugs", "apps_1177_native_frag_indirect_override_unit.wdl")
    // the native app has an instance type of mem1_ssd1_x2, with 2 CPUs and 4 Gb RAM. The tasks request 30 Gb
    val stages = compileGetStages(path.toString, cFlags)
    val stagesSysReq = stages.map { stage =>
      stage.name -> stage.systemRequirements
    }.toMap
    stagesSysReq.size shouldBe 2
    stagesSysReq("if (a)") shouldBe JsObject.empty
    stagesSysReq("apps_1177_mem_int") shouldBe JsObject(
        "*" -> JsObject("instanceType" -> JsString("mem3_ssd1_x4"))
    )
    val stagesExecDetails = stages.map { stage =>
      stage.name -> dxApi.executable(stage.executable).describe(Set(Field.Details))
    }.toMap
    stagesExecDetails("if (a)").details.get.asJsObject.fields
      .get("staticInstanceType")
      .get shouldBe (JsString("mem3_ssd1_x4"))
  }

  it should "add executable links pointing to the intermediate outputs of wrapped workflow" in {
    val path = pathFromBasename("subworkflows", "apps_1175_nested_wf_frag.wdl")
    val stages = compileGetStages(path.toString, cFlags)
    val stagesExecDetails = stages.map { stage =>
      stage.name -> dxApi.executable(stage.executable).describe(Set(Field.Details))
    }.toMap
    val execInfo = stagesExecDetails("frag nested_inner").details.get.asJsObject.fields
      .get("execLinkInfo")
      .get
      .asJsObject
      .fields
      .get("nested_inner")
      .get
      .asJsObject
      .fields
      .get("outputs")
      .getOrElse(JsObject.empty)
    execInfo shouldBe (
        JsObject("test_out" -> JsString("File"),
                 "nested_inner_wf_out" -> JsString("File"))
    )
  }

  it should "compile a workflow with a native app wrapped in frag and keep the default instance" in {
    val path = pathFromBasename("bugs", "apps_1197_native_frag_default_unit.wdl")
    // the native app has an instance type of mem1_ssd1_x2, with 2 CPUs and 4 Gb RAM
    val stages = compileGetStages(path.toString, cFlags)
    val stagesSysReq = stages.map { stage =>
      stage.name -> stage.systemRequirements
    }.toMap
    stagesSysReq.size shouldBe 2
    stagesSysReq("if (a)") shouldBe JsObject.empty
    stagesSysReq("apps_1197_default_instance") shouldBe JsObject.empty
    val stagesExecDetails = stages.map { stage =>
      stage.name -> dxApi.executable(stage.executable).describe(Set(Field.Details))
    }.toMap
    stagesExecDetails("if (a)").details.get.asJsObject.fields
      .get("staticInstanceType")
      .getOrElse(JsObject.empty) shouldBe JsObject.empty
  }

  it should "compile a workflow with a native app and override instance based using RAM and CPU spec" in {
    val path = pathFromBasename("bugs", "apps_1177_native_indirect_override_unit.wdl")
    // the native app has an instance type of mem1_ssd1_x2, with 2 CPUs and 4 Gb RAM. The tasks request 30 Gb and
    // 8 CPU respectively so the instances should be scaled
    val stages = compileGetStages(path.toString, cFlags).map { stage =>
      stage.name -> stage.systemRequirements
    }.toMap
    stages.size shouldBe 3
    stages("default") shouldBe JsObject.empty
    stages("mem_int") shouldBe JsObject(
        "*" -> JsObject("instanceType" -> JsString("mem3_ssd1_x4"))
    )
    stages("cpu_int") shouldBe JsObject(
        "*" -> JsObject("instanceType" -> JsString("mem1_ssd1_x8"))
    )
  }

  it should "compile a workflow with a frag app wrapper using a default instance" in {
    val path = pathFromBasename("frag_runner", "apps_1128_frag_default.wdl")
    val stages = compileGetStages(path.toString, cFlags)
    stages.size shouldBe 1
    val applet = dxApi.applet(stages.filter(_.name.contains("if")).head.executable)
    val instance = applet
      .describe(Set(Field.RunSpec))
      .runSpec
      .get
      .asJsObject
      .fields
      .get("systemRequirements")
      .get
      .asJsObject
      .fields
      .get("main")
      .get
      .asJsObject
      .fields
      .get("instanceType")
    instance.getOrElse("Undefined") shouldBe JsString("mem1_ssd1_v2_x2")
    val sysReq = stages.filter(_.name.contains("if")).head.systemRequirements
    sysReq shouldBe (JsObject.empty)
  }

  it should "compile a task with dynamic instance type selection" in {
    def compile(instanceTypeSelection: String): Map[String, JsValue] = {
      val path = pathFromBasename("compiler", "add.wdl")
      val args = path.toString :: "-instanceTypeSelection" :: instanceTypeSelection :: cFlags
      val retval = Main.compile(args.toVector)
      val appletId = retval match {
        case SuccessfulCompileNativeNoTree(_, Vector(appletId)) => appletId
        case other =>
          throw new Exception(s"expected successful compilation of a single applet not ${other}")
      }
      dxApi.applet(appletId).describe(Set(Field.Details)).details.get.asJsObject.fields
    }
    compile("static") should contain key Constants.InstanceTypeDb
    compile("dynamic") shouldNot contain key Constants.InstanceTypeDb
  }

  it should "Native compile a linear WDL workflow" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "wf_linear.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile a WDL file with a directory input" in {
    val sourceCode = pathFromBasename("input_file", "wdldir.wdl")
    val inputs = pathFromBasename("input_file", "wdldir_input.json")
    val args = List(sourceCode.toString, "--inputs", inputs.toString) ++ cFlags
    val wfId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(id)) => id
      case other                                        => throw new Exception(s"expected success, not ${other}")
    }
    val dxInputsFile = inputs.getParent.resolve(FileUtils.replaceFileSuffix(inputs, ".dx.json"))
    val jsInputs = JsUtils.jsFromFile(dxInputsFile)
    val fields = jsInputs.asJsObject.fields
    fields.size shouldBe 1
    fields("WorkingDir") shouldBe JsString("dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:/tmp/")
    val wf = dxApi.workflow(wfId)
    val params = wf.describe(Set(Field.InputSpec)).inputs.get
    params.size shouldBe 1
    params.head.ioClass shouldBe DxIOClass.String
  }

  it should "Native compile a workflow with a scatter without a call" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "scatter_no_call.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "Native compile a draft2 workflow" taggedAs NativeTest in {
    val path = pathFromBasename("draft2", "shapes.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "handle various conditionals" taggedAs NativeTest in {
    val path = pathFromBasename("draft2", "conditionals_base.wdl")
    val args = path.toString :: cFlags
    /*                :: "--verbose"
            :: "--verboseKey" :: "Native"
            :: "--verboseKey" :: "GenerateIR"*/
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "be able to include pattern information in inputSpec" in {
    val path = pathFromBasename("compiler", "pattern_params.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(Set(Field.InputSpec, Field.OutputSpec))
    val inputParams = desc.inputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    val pattern = inputParams("pattern")
    pattern.help shouldBe Some("The pattern to use to search in_file")
    pattern.group shouldBe Some("Common")
    pattern.label shouldBe Some("Search pattern")
    val inFile = inputParams("in_file")
    inFile.patterns shouldBe Some(IOParameterPatternArray(Vector("*.txt", "*.tsv")))
    inFile.help shouldBe Some("The input file to be searched")
    inFile.group shouldBe Some("Common")
    inFile.label shouldBe Some("Input file")
    val outputParams = desc.outputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    val outFile = outputParams("out_file")
    outFile.patterns shouldBe Some(IOParameterPatternArray(Vector("*.txt", "*.tsv")))
  }

  it should "be able to include pattern object information in inputSpec" in {
    val path = pathFromBasename("compiler", "pattern_obj_params.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(Set(Field.InputSpec, Field.OutputSpec))
    val inputParams = desc.inputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    val pattern = inputParams("pattern")
    pattern.help shouldBe Some("The pattern to use to search in_file")
    pattern.group shouldBe Some("Common")
    pattern.label shouldBe Some("Search pattern")
    val inFile = inputParams("in_file")
    inFile.patterns shouldBe Some(
        IOParameterPatternObject(Some(Vector("*.txt", "*.tsv")),
                                 Some("file"),
                                 Some(Vector("foo", "bar")))
    )
    inFile.help shouldBe Some("The input file to be searched")
    inFile.group shouldBe Some("Common")
    inFile.label shouldBe Some("Input file")
    val outputParams = desc.outputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    val outFile = outputParams("out_file")
    outFile.patterns shouldBe Some(IOParameterPatternArray(Vector("*.txt", "*.tsv")))
  }

  it should "be able to include choices information in inputSpec" in {
    val path = pathFromBasename("compiler", "choice_values.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(Set(Field.InputSpec, Field.OutputSpec))
    val inputParams = desc.inputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    inputParams("pattern").choices shouldBe Some(
        Vector(
            IOParameterValueString(value = "A"),
            IOParameterValueString(value = "B")
        )
    )
    inputParams("in_file").choices shouldBe Some(
        Vector(
            IOParameterValueDataObject("file-Fg5PgBQ0ffP7B8bg3xqB115G"),
            IOParameterValueDataObject("file-Fg5PgBj0ffPP0Jjv3zfv0yxq")
        )
    )
  }

  it should "be able to include annotated choices information in inputSpec" in {
    val path = pathFromBasename("compiler", "choice_obj_values.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val params = inputSpec.inputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    params("pattern").choices shouldBe Some(
        Vector(
            IOParameterValueString(value = "A"),
            IOParameterValueString(value = "B")
        )
    )
    params("in_file").choices shouldBe Some(
        Vector(
            IOParameterValueDataObject(
                id = "file-Fg5PgBQ0ffP7B8bg3xqB115G",
                name = Some("file1")
            ),
            IOParameterValueDataObject(
                id = "file-Fg5PgBj0ffPP0Jjv3zfv0yxq",
                name = Some("file2")
            )
        )
    )
  }

  it should "be able to include suggestion information in inputSpec" in {
    val path = pathFromBasename("compiler", "suggestion_values.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val params = inputSpec.inputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    params("pattern").suggestions shouldBe Some(
        Vector(
            IOParameterValueString("A"),
            IOParameterValueString("B")
        )
    )
    params("in_file").suggestions shouldBe Some(
        Vector(
            IOParameterValueDataObject("file-Fg5PgBQ0ffP7B8bg3xqB115G"),
            IOParameterValueDataObject("file-Fg5PgBj0ffPP0Jjv3zfv0yxq")
        )
    )
  }

  it should "be able to include annotated suggestion information in inputSpec" in {
    val path = pathFromBasename("compiler", "suggestion_obj_values.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val params = inputSpec.inputSpec match {
      case Some(x) => x.map(p => p.name -> p).toMap
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    params("pattern").suggestions shouldBe Some(
        Vector(
            IOParameterValueString("A"),
            IOParameterValueString("B")
        )
    )
    params("in_file").suggestions shouldBe Some(
        Vector(
            IOParameterValueDataObject(
                id = "file-Fg5PgBQ0ffP7B8bg3xqB115G",
                name = Some("file1")
            ),
            DxIoParameterValuePath(
                project = "project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq",
                path = "/test_data/f2.txt.gz",
                name = Some("file2")
            )
        )
    )
  }

  it should "be able to include dx_type information in inputSpec" in {
    val path = pathFromBasename("compiler", "add_dx_type.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val (file_a, file_b) = inputSpec.inputSpec match {
      case Some(x) => (x(0), x(1))
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    file_a.dxType shouldBe Some(DxConstraintString("fastq"))
    file_b.dxType shouldBe Some(
        DxConstraintBool(
            DxConstraintOper.And,
            DxConstraintArray(
                DxConstraintString("fastq"),
                DxConstraintBool(
                    DxConstraintOper.Or,
                    DxConstraintArray(
                        DxConstraintString("Read1"),
                        DxConstraintString("Read2")
                    )
                )
            )
        )
    )
  }

  it should "be able to include default information in inputSpec" in {
    val path = pathFromBasename("compiler", "add_default.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val (int_a, int_b) = inputSpec.inputSpec match {
      case Some(x) => (x(0), x(1))
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    int_a.default shouldBe Some(IOParameterValueNumber(1))
    int_b.default shouldBe Some(IOParameterValueNumber(2))
  }

  it should "be able to include help information in inputSpec" in {
    val path = pathFromBasename("compiler", "add_help.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val (a, b, c, d, e) = inputSpec.inputSpec match {
      case Some(x) => (x(0), x(1), x(2), x(3), x(4))
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    a.help shouldBe Some("lefthand side")
    b.help shouldBe Some("righthand side")
    c.help shouldBe Some("Use this")
    d.help shouldBe Some("Use this")
    e.help shouldBe Some("Use this")
  }

  it should "be able to include group information in inputSpec" in {
    val path = pathFromBasename("compiler", "add_group.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val (a, b) = inputSpec.inputSpec match {
      case Some(x) => (x(0), x(1))
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    a.group shouldBe Some("common")
    b.group shouldBe Some("obscure")
  }

  it should "be able to include label information in inputSpec" in {
    val path = pathFromBasename("compiler", "add_label.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val inputSpec = dxApplet.describe(Set(Field.InputSpec))
    val (a, b) = inputSpec.inputSpec match {
      case Some(x) => (x(0), x(1))
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    a.label shouldBe Some("A positive integer")
    b.label shouldBe Some("A negative integer")
  }

  it should "be able to include information from task meta and extras" in {
    val expectedUpstreamProjects =
      """
        |[
        |  {
        |    "author":"Broad Institute",
        |    "license":"BSD-3-Clause",
        |    "licenseUrl":"https://github.com/broadinstitute/LICENSE.TXT",
        |    "name":"GATK4",
        |    "repoUrl":"https://github.com/broadinstitute/gatk",
        |    "version":"GATK-4.0.1.2"
        |    }
        |]
            """.stripMargin.parseJson

    val expectedWhatsNew =
      """## Changelog
        |### Version 1.1
        |* Added parameter --foo
        |* Added cowsay easter-egg
        |### Version 1.0
        |* Initial version""".stripMargin

    val path = pathFromBasename("compiler", "add_app_meta.wdl")
    val extraPath = pathFromBasename("compiler/extras", "extras_license.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags
    //:: "--verbose"
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(
        Set(
            Field.Description,
            Field.Details,
            Field.DeveloperNotes,
            Field.Properties,
            Field.Summary,
            Field.Tags,
            Field.Title,
            Field.Types
        )
    )

    desc.description.get should startWith(
        "Adds two int together. This app adds together two integers and returns the sum"
    )
    desc.details match {
      case Some(JsObject(fields)) =>
        fields.foreach {
          case ("contactEmail", JsString(value))    => value shouldBe "joe@dev.com"
          case ("upstreamVersion", JsString(value)) => value shouldBe "1.0"
          case ("upstreamAuthor", JsString(value))  => value shouldBe "Joe Developer"
          case ("upstreamUrl", JsString(value))     => value shouldBe "https://dev.com/joe"
          case ("upstreamLicenses", JsArray(array)) => array shouldBe Vector(JsString("MIT"))
          case ("upstreamProjects", array: JsArray) =>
            array shouldBe expectedUpstreamProjects
          case ("whatsNew", JsString(value))                       => value shouldBe expectedWhatsNew
          case (Constants.InstanceTypeDb, JsString(_))             => () // ignore
          case (Constants.Language, JsString(_))                   => () // ignore
          case (Constants.RuntimeAttributes, JsNull | JsObject(_)) => () // ignore
          case (Constants.Version, JsString(_))                    => () // ignore
          case (Constants.Checksum, JsString(_))                   => () // ignore
          case (Constants.SourceCode, JsString(_))                 => () // ignore
          case (Constants.ParseOptions, JsObject(_))               => () // ignore
          case (Constants.DocContents, JsString(_))                => () // ignore
          // old values for sourceCode - can probably delete these
          case ("womSourceCode", JsString(_))                  => () // ignore
          case ("wdlSourceCode", JsString(_))                  => () // ignore
          case (Constants.OriginalName, JsString(_))           => () // ignore
          case (Constants.StaticInstanceType, JsString(value)) => value shouldBe ("mem1_ssd1_x2")
          case other                                           => throw new Exception(s"Unexpected result ${other}")
        }
      case other => throw new Exception(s"Unexpected result ${other}")
    }
    desc.developerNotes shouldBe Some("Check out my sick bash expression! Three dolla signs!!!")
    desc.properties match {
      case Some(m) => m shouldBe Map("baz" -> "blorf")
      case _       => throw new Exception("No properties")
    }
    desc.summary shouldBe Some("Adds two int together")
    desc.tags shouldBe Some(Set("add", "ints", Constants.CompilerTag))
    desc.title shouldBe Some("Add Ints")
    desc.types shouldBe Some(Vector("Adder"))
  }

  it should "be able to include runtime hints" in {
    val path = pathFromBasename("compiler", "add_runtime_hints.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(
        Set(
            Field.Access,
            Field.IgnoreReuse,
            Field.RunSpec
        )
    )

    desc.runSpec match {
      case Some(JsObject(fields)) =>
        fields("executionPolicy") shouldBe JsObject(
            Map(
                "restartOn" -> JsObject(
                    Map(
                        "*" -> JsNumber(1),
                        "UnresponsiveWorker" -> JsNumber(2),
                        "ExecutionError" -> JsNumber(2)
                    )
                ),
                "maxRestarts" -> JsNumber(5)
            )
        )
        fields("timeoutPolicy") shouldBe JsObject(
            Map(
                "*" -> JsObject(
                    Map(
                        "days" -> JsNumber(0),
                        "hours" -> JsNumber(12),
                        "minutes" -> JsNumber(30)
                    )
                )
            )
        )
      case _ => throw new Exception("Missing runSpec")
    }
    desc.access shouldBe Some(
        JsObject(
            Map(
                "network" -> JsArray(Vector(JsString("*"))),
                "developer" -> JsBoolean(true)
            )
        )
    )
    desc.ignoreReuse shouldBe Some(true)
  }

  it should "timeout can be overriden from the extras file" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "add_timeout_override.wdl")
    val extraPath = pathFromBasename("compiler/extras", "short_timeout.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the timeout is what it should be
    val (_, stdout, _) = SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${appId} --json")

    val timeout = stdout.parseJson.asJsObject.fields.get("runSpec") match {
      case Some(JsObject(x)) =>
        x.get("timeoutPolicy") match {
          case None    => throw new Exception("No timeout policy set")
          case Some(s) => s
        }
      case other => throw new Exception(s"Unexpected result ${other}")
    }
    timeout shouldBe JsObject("*" -> JsObject("hours" -> JsNumber(3)))
  }

  it should "be able to include runtime hints and override extras task default" in {
    val path = pathFromBasename("compiler", "add_runtime_hints.wdl")
    val extraPath = pathFromBasename("compiler/extras", "short_timeout.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags

    //:: "--verbose"
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(
        Set(Field.Access, Field.IgnoreReuse, Field.RunSpec)
    )

    desc.runSpec match {
      case Some(JsObject(fields)) =>
        fields("timeoutPolicy") shouldBe JsObject(
            Map(
                "*" -> JsObject(
                    Map(
                        "days" -> JsNumber(0),
                        "hours" -> JsNumber(12),
                        "minutes" -> JsNumber(30)
                    )
                )
            )
        )
      case _ => throw new Exception("Missing runSpec")
    }

    desc.access shouldBe Some(
        JsObject(
            Map(
                "project" -> JsString("ADMINISTER"),
                "allProjects" -> JsString("VIEW"),
                "network" -> JsArray(Vector(JsString("*"))),
                "developer" -> JsBoolean(true)
            )
        )
    )
  }

  it should "be able to include runtime hints with extras per-task override" in {
    val path = pathFromBasename("compiler", "add_runtime_hints.wdl")
    val extraPath = pathFromBasename("compiler/extras", "task_specific_short_timeout.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags
    //:: "--verbose"
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet = dxApi.applet(appId)
    val desc = dxApplet.describe(
        Set(
            Field.RunSpec,
            Field.Access
        )
    )

    // Sometimes the API only returns the fields with non-zero values
    def fillOut(obj: JsValue): JsObject = {
      obj match {
        case JsObject(fields) =>
          val defaults = Map(
              "days" -> JsNumber(0),
              "hours" -> JsNumber(0),
              "minutes" -> JsNumber(0)
          )
          JsObject(fields.view.mapValues {
            case JsObject(inner) => JsObject(defaults ++ inner)
            case _               => throw new Exception("Expected JsObject")
          }.toMap)
        case _ => throw new Exception("Expected JsObject")
      }
    }

    desc.runSpec match {
      case Some(JsObject(fields)) =>
        fillOut(fields("timeoutPolicy")) shouldBe JsObject(
            Map(
                "*" -> JsObject(
                    Map(
                        "days" -> JsNumber(0),
                        "hours" -> JsNumber(3),
                        "minutes" -> JsNumber(0)
                    )
                )
            )
        )
      case _ => throw new Exception("Missing runSpec")
    }
    desc.access shouldBe Some(
        JsObject(
            Map(
                "project" -> JsString("VIEW"),
                "allProjects" -> JsString("VIEW"),
                "network" -> JsArray(Vector(JsString("*"))),
                "developer" -> JsBoolean(true)
            )
        )
    )
  }

  it should "be able to build access for applet with dockerRegistry" in {
    val path = pathFromBasename("compiler", "help_input_params.wdl")
    val extraPath = pathFromBasename("non_spec", "dependency_report_extras.json")
    val args = path.toString :: "-extras" :: extraPath.toString :: cFlags
    //:: "--verbose"
    val (appId1, appId2) = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, x: Vector[String]) => (x.head, x.last)
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxApplet1 = dxApi.applet(appId1)
    val desc1 = dxApplet1.describe(
        Set(
            Field.Access
        )
    )

    desc1.access shouldBe Some(
        JsObject(
            Map(
                "allProjects" -> JsString("VIEW"),
                "network" -> JsArray(Vector.empty)
            )
        )
    )

    val dxApplet2 = dxApi.applet(appId2)
    val desc2 = dxApplet2.describe(
        Set(
            Field.Access
        )
    )

    desc2.access shouldBe Some(
        JsObject(
            Map(
                "allProjects" -> JsString("VIEW"),
                "network" -> JsArray(Vector(JsString("*")))
            )
        )
    )

  }

  it should "be able to build access for tasks in a workflow with dockerRegistry" in {
    val path = pathFromBasename("non_spec", "dependency_report_wf2.wdl")
    val extraPath = pathFromBasename("non_spec", "dependency_report_extras.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags
    //:: "--verbose"
    val wfId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val wf = dxApi.workflow(wfId)
    val wfDesc = wf.describe(Set(Field.Stages))
    val wfStages = wfDesc.stages.get

    val fragAppletId =
      wfStages.collectFirst({ case s if s.name == "frag dependency_report_t1" => s.executable }).get
    val dxFragApplet = dxApi.applet(fragAppletId)
    val fragDesc = dxFragApplet.describe(
        Set(
            Field.Access,
            Field.Details
        )
    )
    fragDesc.access shouldBe Some(
        JsObject(
            Map(
                "allProjects" -> JsString("VIEW"),
                "network" -> JsArray(Vector(JsString("*")))
            )
        )
    )
    val taskAppletId = fragDesc.details match {
      case Some(d) => {
        val execLinkInfo = d.asJsObject.fields.getOrElse(
            Constants.ExecLinkInfo,
            throw new Exception(s"Missing ExecLinkInfo")
        )
        execLinkInfo match {
          case JsObject(e) => {
            val task = e.getOrElse(
                "dependency_report_t1",
                throw new Exception(s"Missing task dependency_report_t1")
            )
            task match {
              case JsObject(t) => {
                val taskId = t.getOrElse(
                    "id",
                    throw new Exception(s"Missing dependency_report_t1 applet id")
                )
                JsUtils.getString(taskId)
              }
              case _ => throw new Exception("Expected dependency_report_t1 applet info")
            }
          }
          case _ => throw new Exception("Expected ExecLinkInfo as a map")
        }
      }
      case _ => throw new Exception(s"Missing details")
    }

    val dxTaskApplet = dxApi.applet(taskAppletId)
    val taskDesc = dxTaskApplet.describe(
        Set(
            Field.Access
        )
    )
    taskDesc.access shouldBe Some(
        JsObject(
            Map(
                "allProjects" -> JsString("VIEW"),
                "network" -> JsArray(Vector(JsString("*")))
            )
        )
    )

    val dxNativeApp = dxApi.resolveApp(
        wfStages.collectFirst({ case s if s.name == "dependency_report_t2" => s.executable }).get
    )
    val nativeDesc = dxNativeApp.describe(
        Set(
            Field.Access
        )
    )
    nativeDesc.access shouldBe Some(
        JsObject(
            Map("network" -> JsArray(Vector.empty))
        )
    )
  }

  it should "be able to include information from workflow meta" in {
    val path = pathFromBasename("compiler", "wf_meta.wdl")
    val args = path.toString :: cFlags
    val wfId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxWorkflow = dxApi.workflow(wfId)
    val desc = dxWorkflow.describe(
        Set(
            Field.Description,
            Field.Details,
            Field.Properties,
            Field.Summary,
            Field.Tags,
            Field.Title,
            Field.Types
        )
    )
    desc.description.get should startWith("This is a workflow that defines some metadata")
    desc.details match {
      case Some(JsObject(fields)) =>
        fields.foreach {
          case ("whatsNew", JsString(value))               => value shouldBe "v1.0: First release"
          case ("delayWorkspaceDestruction", JsBoolean(_)) => ()
          case ("link_inc", JsObject(_))                   => ()
          case ("link_mul", JsObject(_))                   => ()
          case ("link_add", JsObject(_))                   => () // New link added from dependencies. APPS-1175
          case ("execTree", JsString(_))                   => ()
          case (Constants.Version, JsString(_))            => () // ignore
          case (Constants.Checksum, JsString(_))           => () // ignore
          case (Constants.DocContents, JsString(_))        => () // ignore
          case (Constants.SourceCode, JsString(_))         => () // ignore
          case (Constants.ParseOptions, JsObject(_))       => () // ignore
          // old values for sourceCode - can probalby delete these
          case ("womSourceCode", JsString(_))                   => ()
          case ("wdlSourceCode", JsString(_))                   => ()
          case ("staticInstanceTypeSelection", JsBoolean(true)) => ()
          case (Constants.OriginalName, JsString(_))            => ()
          case other                                            => throw new Exception(s"Unexpected result ${other}")
        }
      case other => throw new Exception(s"Unexpected result ${other}")
    }
    desc.properties match {
      case Some(m) => m shouldBe Map("foo" -> "bar")
      case _       => throw new Exception("No properties")
    }
    desc.summary shouldBe Some("A workflow that defines some metadata")
    desc.tags shouldBe Some(Set("foo", "bar", Constants.CompilerTag))
    desc.title shouldBe Some("Workflow with metadata")
    desc.types shouldBe Some(Vector("calculator"))
  }

  it should "be able to include information from workflow parameter meta" in {
    val path = pathFromBasename("compiler", "wf_param_meta.wdl")
    val args = path.toString :: cFlags
    val wfId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val dxWorkflow = dxApi.workflow(wfId)
    val desc = dxWorkflow.describe(Set(Field.Inputs))
    val (x, y) = desc.inputs match {
      case Some(s) => (s(0), s(1))
      case other   => throw new Exception(s"Unexpected result ${other}")
    }
    x.label shouldBe Some("Left-hand side")
    x.default shouldBe Some(IOParameterValueNumber(3))
    y.label shouldBe Some("Right-hand side")
    y.default shouldBe Some(IOParameterValueNumber(5))
  }

  it should "dependency report should be added to WF description" in {
    val path = pathFromBasename("non_spec", "dependency_report_wf1.wdl")
    val extraPath = pathFromBasename("non_spec", "dependency_report_extras.json")
    val args = path.toString :: "-extras" :: extraPath.toString :: cFlags
    val wfId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"Unexpected compilation result ${other}")
    }

    val dxWorkflow = dxApi.workflow(wfId)
    val desc = dxWorkflow.describe(Set(Field.Description))
    desc.description match {
      case Some(d) => {
        d should include("app-BZ9ZQzQ02VP5gpkP54b96pYY")
        d should include("file-G4BV6180yzZyvZ124KB0q46P")
        d should include("file-G4BV61j0yzZgf6JQKxP4gQ3Y")
        d should include("file-G4BV61Q0yzZq60Jj4K5vfG92")
        d should include("file-G4BV6280yzZq60Jj4K5vfG96")
        d should include("file-G4BV6200yzZq60Jj4K5vfG94")
        d should include("file-G4BV6280yzZkpgjj4Jx6Fjj9")
        d should include("file-G4BV6100yzZz17bx4JkQkybb")
        d should include("alpine:3.14")
        d should include("ubuntu:20.04")
        d should include("mem1_ssd2_x4")
        d should include(
            "Available instance types determined during WDL compilation will be used"
        )
      }
      case _ => throw new Exception(s"Expected ${wfId} to have a description.")
    }
  }

  private def validateBundledDepends(appletId: String, bundledName: String): Unit = {
    val applet = dxApi.applet(appletId)
    val desc = applet.describe(Set(Field.RunSpec))

    desc.runSpec match {
      case Some(rs) => {
        val bundledDepends = rs.asJsObject.fields.getOrElse(
            "bundledDepends",
            throw new Exception(s"Expected ${appletId} to have bundledDepends")
        )
        bundledDepends match {
          case JsArray(items) => {
            val searchItem = items.find(v => v.asJsObject.toString.contains(bundledName))
            searchItem should not be empty
          }
          case _ => throw new Exception("Expected bundledDepends to be array")
        }
      }
      case _ => throw new Exception(s"Expected ${appletId} to have runSpec")
    }
  }

  it should "add clonable dependencies to bundledDepends" in {
    val path = pathFromBasename("non_spec", "apps_623_wf.wdl")
    val args = path.toString :: cFlags
    val wfId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"Unexpected compilation result ${other}")
    }

    // Describe workflow, get stages
    val wf = dxApi.workflow(wfId)
    val desc = wf.describe(Set(Field.Stages))

    // Validate bundledDepends contents of applets
    desc.stages match {
      case Some(s) => {
        val t1_frag = s.find(stageDesc => stageDesc.id == "stage-1")
        validateBundledDepends(
            t1_frag
              .getOrElse(throw new Exception(s"Expected to find stage-1 in ${wfId}"))
              .executable,
            "apps_623_t1"
        )
        val t2_frag = s.find(stageDesc => stageDesc.id == "stage-3")
        validateBundledDepends(
            t2_frag
              .getOrElse(throw new Exception(s"Expected to find stage-3 in ${wfId}"))
              .executable,
            "apps_623_t2"
        )
        val t3 = s.find(stageDesc => stageDesc.id == "stage-4")
        validateBundledDepends(
            t3.getOrElse(throw new Exception(s"Expected to find stage-4 in ${wfId}")).executable,
            "ubuntu_20_04.tar.gz"
        )
      }
      case _ => throw new Exception(s"Expected ${wfId} to have stages")
    }
  }

  it should "deep nesting" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "environment_passing_deep_nesting.wdl")
    val args = path.toString :: cFlags
    /*                :: "--verbose"
            :: "--verboseKey" :: "Native"
            :: "--verboseKey" :: "GenerateIR"*/
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "make default task timeout 48 hours" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "add_timeout.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the timeout is what it should be
    val (_, stdout, _) = SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${appId} --json")

    val timeout = stdout.parseJson.asJsObject.fields.get("runSpec") match {
      case Some(JsObject(x)) =>
        x.get("timeoutPolicy") match {
          case None    => throw new Exception("No timeout policy set")
          case Some(s) => s
        }
      case other => throw new Exception(s"Unexpected result ${other}")
    }
    timeout shouldBe JsObject(
        "*" -> JsObject("days" -> JsNumber(2), "hours" -> JsNumber(0), "minutes" -> JsNumber(0))
    )
  }

  it should "allow choosing GPU instances" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "GPU2.wdl")
    val args = path.toString :: cFlags
    val appId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the timeout is what it should be
    val (_, stdout, _) = SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${appId} --json")
    val obj = stdout.parseJson.asJsObject
    val obj2 = obj.fields("runSpec").asJsObject
    val obj3 = obj2.fields("systemRequirements").asJsObject
    val obj4 = obj3.fields("main").asJsObject
    val instanceType = obj4.fields.get("instanceType") match {
      case Some(JsString(x)) => x
      case other             => throw new Exception(s"Unexpected result ${other}")
    }
    instanceType should include("_gpu")
  }

  it should "Compile a workflow with subworkflows on the platform with the reorg app" in {
    val path = pathFromBasename("subworkflows", basename = "trains_station.wdl")
    val extrasContent =
      s"""|{
          | "custom_reorg" : {
          |    "app_id" : "${reorgAppletId.get}",
          |    "conf" : null
          |  }
          |}
          |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "-extras" :: extrasPath :: cFlagsReorgCompile
      Main.compile(args.toVector)
    }
    val wfId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val wf = dxApi.workflow(wfId)
    val wfDesc = wf.describe(Set(Field.Stages))
    val wfStages = wfDesc.stages.get

    // there should be 4 stages: 1) common 2) train_stations 3) outputs 4) reorg
    wfStages.size shouldBe 4
    val reorgStage = wfStages.last

    reorgStage.id shouldBe "stage-reorg"
    reorgStage.executable shouldBe reorgAppletId.get

    // There should be 3 inputs, the output from output stage and the custom reorg config file.
    val reorgInput: JsObject = reorgStage.input match {
      case JsObject(x) => JsObject(x)
      case _           => throw new Exception("unexpected")
    }

    // no reorg conf input. only status.
    reorgInput.fields.size shouldBe 1
    reorgInput.fields.keys should contain theSameElementsAs Set(Constants.ReorgStatus.encoded)
  }

  // ignore for now as the test will fail in staging
  it should "Compile a workflow with subworkflows on the platform with the reorg app with config file in the input" in {
    // This works in conjunction with "Compile a workflow with sub-workflows on the platform with the reorg app".
    val path = pathFromBasename("subworkflows", basename = "trains_station.wdl")
    // upload random file
    val (_, uploadOut, _) = SysUtils.execCommand(
        s"dx upload ${path.toString} --destination /reorg_tests --brief"
    )
    val fileId = uploadOut.trim
    val extrasContent =
      s"""|{
          | "custom_reorg" : {
          |    "app_id" : "${reorgAppletId.get}",
          |    "conf" : "dx://$fileId"
          |  }
          |}
          |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "-extras" :: extrasPath :: cFlagsReorgCompile
      Main.compile(args.toVector)
    }
    val wfId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    val wf = dxApi.workflow(wfId)
    val wfDesc = wf.describe(Set(Field.Stages))
    val wfStages = wfDesc.stages.get
    val reorgStage = wfStages.last

    // There should be 3 inputs, the output from output stage and the custom reorg config file.
    val reorgInput: JsObject = reorgStage.input match {
      case JsObject(x) => JsObject(x)
      case _           => throw new Exception("unexpected")
    }
    // no reorg conf input. only status.
    reorgInput.fields.size shouldBe 2
    reorgInput.fields.keys should contain theSameElementsAs Set(Constants.ReorgStatus.encoded,
                                                                Constants.ReorgConfig.encoded)
  }

  it should "ensure subworkflow with custom reorg app does not contain reorg attribute" in {
    val path = pathFromBasename("subworkflows", basename = "trains_station.wdl")
    // upload random file
    val extrasContent =
      s"""|{
          | "custom_reorg" : {
          |    "app_id" : "${reorgAppletId.get}",
          |    "conf" : null
          |  }
          |}
          |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "-extras" :: extrasPath :: cFlagsReorgIR
      Main.compile(args.toVector)
    }
    val bundle = retval match {
      case SuccessfulCompileIR(bundle) => bundle
      case _                           => throw new Exception("unexpected")
    }

    // this is a subworkflow so there is no reorg_status___ added.
    val trainsOutputVector: Callable = bundle.allCallables("trains")
    trainsOutputVector.outputVars.size shouldBe 1
  }

  it should "compile workflow with task imported by multiple paths" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", basename = "multi_import.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile workflow with call to task having optional output" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", basename = "optional_call_output.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile workflow with optional output" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", basename = "optional_output.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "set job-reuse flag on applet" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "add_timeout.wdl")
    val extrasContent =
      """|{
         |  "ignoreReuse": true
         |}
         |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "--extras" :: extrasPath :: cFlags
      Main.compile(args.toVector)
    }
    val appletId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the job reuse flag is set
    val (_, stdout, _) =
      SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${appletId} --json")
    val ignoreReuseFlag = stdout.parseJson.asJsObject.fields.get("ignoreReuse")
    ignoreReuseFlag shouldBe Some(JsBoolean(true))
  }

  it should "override job-reuse flag on applet" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "add_runtime_hints.wdl")
    val extrasContent =
      """|{
         |  "ignoreReuse": false
         |}
         |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "--extras" :: extrasPath :: cFlags
      Main.compile(args.toVector)
    }
    val appletId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the job reuse flag is set
    val (_, stdout, _) =
      SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${appletId} --json")
    val ignoreReuseFlag = stdout.parseJson.asJsObject.fields.get("ignoreReuse")
    ignoreReuseFlag shouldBe Some(JsBoolean(false))
  }

  it should "set job-reuse flag on workflow" taggedAs NativeTest in {
    val path = pathFromBasename("subworkflows", basename = "check_route.wdl")
    val extrasContent =
      """|{
         |  "ignoreReuse": true
         |}
         |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "-extras" :: extrasPath :: cFlags
      Main.compile(args.toVector)
    }
    val wfId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the job reuse flag is set
    val (_, stdout, _) =
      SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${wfId} --json")
    val wfDesc = stdout.parseJson.asJsObject
    val ignoreReuseFlag = wfDesc.fields.get("ignoreReuse")
    ignoreReuseFlag shouldBe Some(JsArray(JsString("*")))

    val taskIgnoreReuseFlag = {
      val taskId = wfDesc.fields
        .get("links")
        .map(
            {
              case l: JsArray => l.elements.map(_.toString).head
              case _          => None
            }
        )
        .getOrElse("Concat")
      val (_, stdout, _) =
        SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${taskId} --json")
      stdout.parseJson.asJsObject.fields.get("ignoreReuse")
    }
    taskIgnoreReuseFlag shouldBe Some(JsBoolean(true))
  }

  it should "set delayWorkspaceDestruction on applet" taggedAs NativeTest in {
    val path = pathFromBasename("compiler", "add_timeout.wdl")
    val extrasContent =
      """|{
         |  "delayWorkspaceDestruction": true
         |}
         |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "-extras" :: extrasPath :: cFlags
      Main.compile(args.toVector)
    }
    val appletId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the delayWorkspaceDestruction flag is set
    val (_, stdout, _) =
      SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${appletId} --json")
    val details = stdout.parseJson.asJsObject.fields("details")
    val delayWD = details.asJsObject.fields.get("delayWorkspaceDestruction")
    delayWD shouldBe Some(JsTrue)
  }

  it should "set delayWorkspaceDestruction on workflow" taggedAs NativeTest in {
    val path = pathFromBasename("subworkflows", basename = "trains_station.wdl")
    val extrasContent =
      """|{
         |  "delayWorkspaceDestruction": true
         |}
         |""".stripMargin
    val retval = WithExtras(extrasContent) { extrasPath =>
      val args = path.toString :: "-extras" :: extrasPath :: cFlags
      Main.compile(args.toVector)
    }
    val wfId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }

    // make sure the flag is set on the resulting workflow
    val (_, stdout, _) =
      SysUtils.execCommand(s"dx describe ${dxTestProject.id}:${wfId} --json")
    val details = stdout.parseJson.asJsObject.fields("details")
    val delayWD = details.asJsObject.fields.get("delayWorkspaceDestruction")
    delayWD shouldBe Some(JsTrue)
  }

  it should "Native compile a CWL tool" taggedAs NativeTest in {
    val path = pathFromBasename("cwl", "cat.cwl.json")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "Compile a tool with -useManifests flag" in {
    val path = pathFromBasename("compiler", "add.wdl")
    val args = path.toString :: "-useManifests" :: cFlags
    val appletId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(x)) => x
      case other =>
        throw new AssertionError(s"expected Success, got ${other}")
    }
    val dxApplet = dxApi.applet(appletId)
    val desc = dxApplet.describe()

    val input = desc.inputSpec.get.map(i => i.name -> i).toMap
    input.keySet shouldBe Set(
        Constants.InputManifestFiles,
        Constants.InputManifest,
        Constants.InputManifest.addSuffix(Constants.FlatFilesSuffix),
        Constants.InputLinks,
        Constants.InputLinks.addSuffix(Constants.FlatFilesSuffix),
        Constants.WorkflowInputManifest,
        Constants.WorkflowInputManifest.addSuffix(Constants.FlatFilesSuffix),
        Constants.WorkflowInputManifestFiles,
        Constants.WorkflowInputLinks,
        Constants.WorkflowInputLinks.addSuffix(Constants.FlatFilesSuffix),
        Constants.OutputId,
        Constants.ExtraOutputs,
        Constants.ExtraOutputs.addSuffix(Constants.FlatFilesSuffix),
        Constants.CallName,
        Constants.Overrides,
        Constants.Overrides.addSuffix(Constants.FlatFilesSuffix)
    ).map(_.encoded)
    input(Constants.InputManifest.encoded).ioClass shouldBe DxIOClass.Hash
    input(Constants.InputManifest.encoded).optional shouldBe true
    input(Constants.InputManifestFiles.encoded).ioClass shouldBe DxIOClass.FileArray
    input(Constants.InputLinks.encoded).ioClass shouldBe DxIOClass.Hash
    input(Constants.InputLinks.encoded).optional shouldBe true
    input(Constants.WorkflowInputManifest.encoded).ioClass shouldBe DxIOClass.Hash
    input(Constants.WorkflowInputManifest.encoded).optional shouldBe true
    input(Constants.WorkflowInputManifestFiles.encoded).ioClass shouldBe DxIOClass.FileArray
    input(Constants.WorkflowInputLinks.encoded).ioClass shouldBe DxIOClass.Hash
    input(Constants.WorkflowInputLinks.encoded).optional shouldBe true
    input(Constants.OutputId.encoded).ioClass shouldBe DxIOClass.String
    input(Constants.ExtraOutputs.encoded).ioClass shouldBe DxIOClass.Hash
    input(Constants.ExtraOutputs.encoded).optional shouldBe true
    input(Constants.CallName.encoded).ioClass shouldBe DxIOClass.String
    input(Constants.CallName.encoded).optional shouldBe true

    val output = desc.outputSpec.get.map(i => i.name -> i).toMap
    output.keySet shouldBe Set(Constants.OutputManifest.encoded)
    output(Constants.OutputManifest.encoded).ioClass shouldBe DxIOClass.File
  }

  it should "compile a workflow with a native app with file output" taggedAs NativeTest in {
    val path = pathFromBasename("bugs", "native_with_file_output.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    val wfId = retval match {
      case SuccessfulCompileNativeNoTree(_, Vector(wfId)) => wfId
      case other                                          => throw new Exception(s"expected single workflow not ${other}")
    }
    // the native app has an instance type of x2, but the WDL task specifies dx_instance_type
    // of x8, so make sure the stage overrides the instance type
    val stages = dxApi.workflow(wfId).describe().stages.get
    stages.size shouldBe 2
    stages.head.executable shouldBe "app-aws_s3_to_platform_files/1.0.2"
    stages.head.systemRequirements shouldBe JsObject(
        "*" -> JsObject("instanceType" -> JsString("mem1_ssd1_v2_x8"))
    )
  }

  it should "compile a packed CWL workflow" in {
    val path = pathFromBasename("cwl", "basename-fields-test.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) match {
      case _: SuccessfulCompileNativeNoTree => ()
      case other =>
        throw new Exception(s"expected success not ${other}")
    }
  }

  it should "compile a workflow with only an output section" in {
    val path = pathFromBasename("draft2", "output_only.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile a CWL workflow with a file default" in {
    val path = pathFromBasename("cwl", "count-lines9-wf-noET.cwl.json")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile a CWL workflow with a stage input unconnected to process input" in {
    val path = pathFromBasename("cwl", "pass-unconnected.cwl.json")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile a CWL workflow with downcasted tool output" in {
    val path = pathFromBasename("cwl", "count-lines11-null-step-wf-noET.cwl.json")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "compile a CWL workflow with a scatter" in {
    val path = pathFromBasename("cwl", "scatter-wf4.cwl.json")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileNativeNoTree]
  }

  it should "archive an identical task" in {
    val folder = s"${unitTestsPath}/testArchive"
    val flags = cFlagsBase ++ List("-compileMode",
                                   "NativeWithoutRuntimeAsset",
                                   "-folder",
                                   folder,
                                   "-locked")
    val path1 = pathFromBasename("compiler", "add.wdl")
    val args1 = path1.toString :: flags
    // compile once
    val appletId1 = Main.compile(args1.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(appletId)) => appletId
      case other =>
        throw new Exception(s"expected single applet not ${other}")
    }
    // compile slightly changed version with archive flag
    val path2 = pathFromBasename("compiler", "add_changed.wdl")
    val args2 = path2.toString :: "-archive" :: flags
    val appletId2 = Main.compile(args2.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(appletId)) => appletId
      case other =>
        throw new Exception(s"expected single applet not ${other}")
    }
    appletId1 should not be appletId2
    // check that first applet has moved to a .archive folder and tagged with "archived"
    val desc = dxApi.applet(appletId1).describe(Set(Field.Tags))
    desc.folder should startWith(s"${folder}/.archive")
    desc.tags should not be empty
    desc.tags.get should contain(DxExecutableDirectory.ArchivedTag)
  }

  it should "reuse identical tasks" in {
    val path = pathFromBasename("bugs", "apps-788.wdl")
    val args = path.toString :: cFlags
    val appletId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(appletId)) => appletId
      case other =>
        throw new Exception(s"expected single applet not ${other}")
    }
    // compiling a second time into a different folder with -projectWideReuse should reuse the same applet
    val args2 = path.toString :: cFlagsReuse
    val appletId2 = Main.compile(args2.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(appletId)) => appletId
      case other =>
        throw new Exception(s"expected single applet not ${other}")
    }
    appletId shouldBe appletId2
  }

  it should "reuse identical cwl tasks" in {
    val path = pathFromBasename("cwl", "params.cwl.json")
    val args = path.toString :: cFlags
    val appletId = Main.compile(args.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(appletId)) => appletId
      case other =>
        throw new Exception(s"expected single applet not ${other}")
    }
    // compiling a second time into a different folder with -projectWideReuse should reuse the same applet
    val args2 = path.toString :: cFlagsReuse
    val appletId2 = Main.compile(args2.toVector) match {
      case SuccessfulCompileNativeNoTree(_, Vector(applet2Id)) => applet2Id
      case other =>
        throw new Exception(s"expected single applet not ${other}")
    }
    appletId shouldBe appletId2
  }
//  it should "compile a task with a string + int concatenation" taggedAs NativeTest in {
//    val path = pathFromBasename("non_spec", "string_int_concat.wdl")
//    val args = path.toString :: "-wdlMode" :: "lenient" :: cFlags
//    val retval = Main.compile(args.toVector)
//    retval shouldBe a[SuccessfulCompileNativeNoTree]
//  }
}
