package dx.dxni

import java.nio.file.{Files, Path, Paths}
import dx.Assumptions.{isLoggedIn, toolkitCallable}
import dx.Tags.NativeTest
import dx.api.{DxApi, DxApplet, DxPath}
import dx.core.languages.wdl.WdlUtils
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import dx.util.{FileUtils, Logger, SysUtils}
import dxCompiler.Main
import dxCompiler.Main.{SuccessfulCompileNativeNoTree, SuccessfulDxNI}
import wdlTools.syntax.NoSuchParserException

import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import java.util.UUID.randomUUID

class DxNativeInterfaceTest extends AnyFlatSpec with Matchers with BeforeAndAfterAll {
  assume(isLoggedIn)
  assume(toolkitCallable)
  private val logger = Logger.Quiet
  private val dxApi = DxApi()(logger)

  private val testProject = "dxCompiler_playground"
  private val dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm")
  private val test_time = dateFormatter.format(LocalDateTime.now)

  private lazy val dxTestProject =
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
  private val unitTestsPath = s"unit_tests/${username}"
  private val folderPath =
    s"/${unitTestsPath}/applets_${test_time}_${randomUUID().toString.substring(24)}/"
  private val compilePath = s"${folderPath}compiled/"
  private val logPath: Path = Files.createTempFile("dxni", ".log")

  override def beforeAll(): Unit = {
    // build the directory with the native applets
    dxTestProject.newFolder(folderPath, parents = true)
    // building necessary applets before starting the tests
    val nativeApplets = Vector(
        "native_concat",
        "native_diff",
        "native_mk_list",
        "native_sum",
        "native_sum_012"
    )
    val topDir = Paths.get(System.getProperty("user.dir"))
    nativeApplets.foreach { app =>
      SysUtils.execCommand(
          s"dx build $topDir/test/applets/$app --destination ${testProject}:${folderPath}"
      )
    }
  }

  override def afterAll(): Unit = {
    dxTestProject.removeFolder(folderPath, recurse = true)
    if (Files.exists(logPath)) {
      Files.delete(logPath)
    }
  }

  private def parseWdlTasks(
      wfSource: String
  ): (Map[String, TAT.Task], Map[String, WdlTypes.T], TAT.Document) = {
    val (tDoc, typeAliases) = WdlUtils.parseAndCheckSourceString(wfSource, "test")
    val tasks = tDoc.elements.collect {
      case task: TAT.Task => task.name -> task
    }.toMap
    (tasks, typeAliases.toMap, tDoc)
  }

  private def runDxni(args: Vector[String]): Map[String, TAT.Task] = {
    val outputPath: Path = Files.createTempFile("dx_extern_one", ".wdl")
    try {
      Main.dxni(args ++ Vector("-output", outputPath.toString)) shouldBe a[SuccessfulDxNI]
      // check that the generated file contains the correct tasks
      val content = FileUtils.readFileContent(outputPath)
      val (tasks, _, _) = parseWdlTasks(content)
      tasks
    } catch {
      case _: NoSuchParserException => Map.empty
    } finally {
      if (Files.exists(outputPath)) {
        Files.delete(outputPath)
      }
    }
  }

  private def runCompile(args: Vector[String]): String = {
    val wfId = Main.compile(args) match {
      case SuccessfulCompileNativeNoTree(_, Vector(wfId)) => wfId
      case other                                          => throw new Exception(s"expected single workflow not ${other}")
    }
    wfId
  }

  private def pathFromBasename(dir: String, basename: String): Path = {
    val p = getClass.getResource(s"/${dir}/${basename}").getPath
    Paths.get(p)
  }

  it should "be able to build interfaces to native applets" taggedAs NativeTest in {
    val args = Vector("-force",
                      "-quiet",
                      "-folder",
                      folderPath,
                      "-project",
                      dxTestProject.id,
                      "-language",
                      "wdl_1_0",
                      "-apps",
                      "exclude")
    val tasks = runDxni(args)
    tasks.keySet shouldBe Set(
        "native_sum",
        "native_sum_012",
        "native_mk_list",
        "native_diff",
        "native_concat"
    )
  }

  it should "be able to build an interface to a specific applet" taggedAs NativeTest in {
    val args = Vector(
        "-force",
        "-quiet",
        "-path",
        s"${folderPath}/native_sum",
        "-project",
        dxTestProject.id,
        "-language",
        "wdl_1_0",
        "-apps",
        "exclude"
    )
    val tasks = runDxni(args)
    tasks.keySet shouldBe Set("native_sum")
  }

  it should "be able to build an interface to a specific app" taggedAs NativeTest in {
    val args = Vector(
        "-force",
        "-quiet",
        "-path",
        "app-qualimap2_anlys",
        "-language",
        "wdl_1_0"
    )
    val tasks = runDxni(args)
    tasks.keySet shouldBe Set("qualimap2_anlys")
  }

  it should "build an interface to an applet specified by ID" taggedAs NativeTest in {
    val applet = dxApi.resolveDataObject(
        DxPath.format(dxTestProject.id, folderPath, "native_sum")
    ) match {
      case applet: DxApplet => applet
      case other            => throw new Exception(s"${other} not an applet")
    }
    val args = Vector("-force",
                      "-quiet",
                      "-path",
                      applet.id,
                      "-project",
                      dxTestProject.id,
                      "-language",
                      "wdl_1_0",
                      "-apps",
                      "exclude")
    val tasks = runDxni(args)
    tasks.keySet shouldBe Set("native_sum")
  }

  it should "not build any interfaces for dxCompiler applets and throw a warning" taggedAs NativeTest in {
    val path = pathFromBasename("bugs", "apps_1351.wdl").toString
    val compileArgs = Vector(
        "-project",
        dxTestProject.id,
        "-force",
        "-compileMode",
        "NativeWithoutRuntimeAsset",
        "-folder",
        compilePath
    )
    val _ = runCompile(path +: compileArgs)
    val dxniArgs = Vector(
        "-force",
        "-folder",
        compilePath,
        "-project",
        dxTestProject.id,
        "-apps",
        "exclude",
        "-r",
        "-logFile",
        logPath.toString
    )
    val tasks = runDxni(dxniArgs)
    tasks.keySet shouldBe Set.empty
    Files.readString(logPath) should include(
        "Applets apps_1351_outputs, apps_1351_common, aa were ignored by dxni because they were compiled with dxCompiler"
    )
  }
}
