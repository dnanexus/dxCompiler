package dx.executor.wdl

import java.nio.file.{Files, Path, Paths}
import Assumptions.{dxdaCallable, isLoggedIn}
import Tags.{ApiTest, EdgeTest}
import dx.api.{
  DiskType,
  DxAnalysis,
  DxApi,
  DxFile,
  DxFileDescCache,
  DxInstanceType,
  DxJob,
  DxPath,
  DxProject,
  ExecutionEnvironment,
  FileUpload,
  InstanceTypeDB
}
import dx.core.Constants
import dx.core.io.{DxWorkerPaths, StreamFiles}
import dx.core.ir.{
  DxName,
  Manifest,
  ParameterLinkDeserializer,
  ParameterLinkSerializer,
  Type,
  Value
}
import dx.core.languages.wdl.{CodeGenerator, VersionSupport, WdlDxName, WdlOptions, WdlUtils}
import dx.executor.{JobMeta, TaskExecutor, TaskExecutorResult}
import dx.util.{CodecUtils, FileSourceResolver, FileUtils, JsUtils, Logger, PosixPath, SysUtils}
import dx.util.protocols.DxFileAccessProtocol
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._
import wdlTools.types.{TypedAbstractSyntax => TAT}

import scala.util.Random

private case class TaskTestJobMeta(override val workerPaths: DxWorkerPaths,
                                   override val dxApi: DxApi = DxApi.get,
                                   override val logger: Logger = Logger.get,
                                   override val rawJsInputs: Map[DxName, JsValue],
                                   rawInstanceTypeDb: InstanceTypeDB,
                                   rawSourceCode: String,
                                   pathToDxFile: Map[Path, DxFile] = Map.empty,
                                   useManifestInputs: Boolean = false,
                                   downloadFiles: Boolean = true)
    extends JobMeta(workerPaths, WdlDxName, dxApi, logger) {
  var outputs: Option[Map[DxName, JsValue]] = None

  override val project: DxProject = null

  override def uploadFiles(filesToUpload: Iterable[FileUpload]): Map[Path, DxFile] = {
    filesToUpload.map { upload =>
      pathToDxFile
        .collectFirst {
          case (path, dxFile) if upload.source.endsWith(path) => upload.source -> dxFile
        }
        .getOrElse {
          throw new Exception(
              s"${upload.source} does not match any of ${pathToDxFile.keys.mkString(",")}"
          )
        }
    }.toMap
  }

  override def writeRawJsOutputs(outputJs: Map[DxName, JsValue]): Unit = {
    outputs = Some(outputJs)
  }

  override lazy val jobId: String = s"job-${Random.alphanumeric.take(24).mkString}"

  override def runJobScriptFunction(name: String,
                                    successCodes: Option[Set[Int]] = Some(Set(0)),
                                    truncateLogs: Boolean = false): Unit = {
    name match {
      case TaskExecutor.DownloadDxda if downloadFiles && !dxdaCallable =>
        throw new Exception("cannot call dxda")
      case TaskExecutor.DownloadDxda if downloadFiles =>
        val dxdaManifest = workerPaths.getDxdaManifestFile().toString
        SysUtils.execCommand(
            s"""bzip2 ${dxdaManifest} &&
               |cd / &&
               |dx-download-agent download ${dxdaManifest}.bz2""".stripMargin
              .replaceAll("\n", " ")
        )
      case TaskExecutor.RunCommand =>
        val script: Path = workerPaths.getCommandFile().asJavaPath
        if (Files.exists(script)) {
          val (_, stdout, stderr) =
            SysUtils.execCommand(script.toString, exceptionOnFailure = false)
          val rc = FileUtils.readFileContent(workerPaths.getReturnCodeFile().asJavaPath).trim.toInt
          if (!successCodes.forall(_.contains(rc))) {
            throw new Exception(
                s"job script ${script} failed with return code ${rc}\nstdout:\n${stdout}\nstderr:\n${stderr}"
            )
          }
        }
      case _ => ()
    }
  }

  override val analysis: Option[DxAnalysis] = None

  override val parentJob: Option[DxJob] = None

  override val instanceType: Option[String] = Some(TaskTestJobMeta.InstanceType)

  override def folder: String = "/"

  override def getJobDetail(name: String): Option[JsValue] = None

  override def getExecutableAttribute(name: String): Option[JsValue] = None

  private val executableDetails: Map[String, JsValue] = Map(
      Constants.InstanceTypeDb -> JsString(
          CodecUtils.gzipAndBase64Encode(
              rawInstanceTypeDb.toJson.prettyPrint
          )
      ),
      Constants.SourceCode -> JsString(CodecUtils.gzipAndBase64Encode(rawSourceCode)),
      Constants.UseManifests -> JsBoolean(useManifestInputs)
  )

  override def getExecutableDetail(name: String): Option[JsValue] = {
    executableDetails.get(name)
  }

  override def error(e: Throwable): Unit = {}
}

private object TaskTestJobMeta {
  val InstanceType = "mem_ssd_unicorn"
}

// This test module requires being logged in to the platform.
class WdlTaskExecutorTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)
  private val logger = Logger.Verbose
  private val dxApi = DxApi()(logger)
  private val unicornInstance = DxInstanceType(
      TaskTestJobMeta.InstanceType,
      100,
      100,
      4,
      gpu = false,
      Vector(
          ExecutionEnvironment(Constants.OsDistribution,
                               Constants.OsRelease,
                               Vector(Constants.OsVersion))
      ),
      Some(DiskType.SSD),
      Some(1)
  )
  private val instanceTypeDB =
    InstanceTypeDB(Map(TaskTestJobMeta.InstanceType -> unicornInstance))

  private def pathFromBasename(basename: String): Option[Path] = {
    getClass.getResource(s"/task_runner/${basename}") match {
      case null => None
      case res  => Some(Paths.get(res.getPath))
    }
  }

  // Recursively go into a wdlValue, and add a base path to the file.
  // For example:
  //   foo.txt ---> /home/joe_heller/foo.txt
  //
  // This is used to convert relative paths to test files into absolute paths.
  // For example, convert:
  //  {
  //    "pattern" : "snow",
  //    "in_file" : "manuscript.txt"
  //  }
  // into:
  //  {
  //    "pattern" : "snow",
  //    "in_file" : "/home/joe_heller/dxCompiler/src/test/resources/runner_tasks/manuscript.txt"
  // }
  //
  object AddBaseDir extends Value.TransformHandler {
    override def apply(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      (t, value) match {
        // single file
        case (Some(Type.TFile) | None, f: Value.VFile) if !f.uri.startsWith(DxPath.DxUriPrefix) =>
          pathFromBasename(f.uri).map(path => f.copy(uri = path.toString)).orElse {
            throw new Exception(s"File ${f.uri} does not exist")
          }
        case (Some(Type.TFile) | None, Value.VString(uri)) if !uri.startsWith(DxPath.DxUriPrefix) =>
          pathFromBasename(uri).map(path => Value.VFile(uri = path.toString)).orElse {
            throw new Exception(s"File ${uri} does not exist")
          }
        case _ => None
      }
    }
  }

  private def getInputs(wdlName: String): Map[DxName, JsValue] = {
    pathFromBasename(s"${wdlName}_input.json") match {
      case Some(path) if Files.exists(path) =>
        JsUtils.getFields(JsUtils.jsFromFile(path)).map {
          case (name, jsv) => WdlDxName.fromEncodedName(name) -> jsv
        }
      case _ => Map.empty
    }
  }

  private def getExpectedOutputs(wdlName: String): Option[Map[DxName, JsValue]] = {
    pathFromBasename(s"${wdlName}_output.json") match {
      case Some(path) if Files.exists(path) =>
        Some(JsUtils.getFields(JsUtils.jsFromFile(path)).map {
          case (name, jsv) => WdlDxName.fromDecodedName(name) -> jsv
        })
      case _ => None
    }
  }

  private def createTaskExecutor(
      wdlName: String,
      streamFiles: StreamFiles.StreamFiles = StreamFiles.None,
      pathToDxFile: Map[Path, DxFile] = Map.empty,
      useManifests: Boolean = false,
      downloadFiles: Boolean = true
  ): (WdlTaskExecutor, TaskTestJobMeta) = {
    val wdlFile: Path = pathFromBasename(s"${wdlName}.wdl").get
    val inputs = getInputs(wdlName)
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()

    // create FileSourceResolver
    val dxFileDescCache = DxFileDescCache.empty
    val inputDeserializer: ParameterLinkDeserializer =
      ParameterLinkDeserializer(dxFileDescCache, dxApi)
    val dxProtocol = DxFileAccessProtocol(dxApi, dxFileDescCache)
    val fileResolver = FileSourceResolver.create(
        localDirectories = Vector(jobRootDir),
        userProtocols = Vector(dxProtocol),
        logger = logger
    )

    // create a stand-alone task
    val wdlOptions = WdlOptions.default
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceFile(wdlFile, wdlOptions, fileResolver)
    val codegen = CodeGenerator(typeAliases.toMap, doc.version.value, logger)
    val tasks: Vector[TAT.Task] = doc.elements.collect {
      case task: TAT.Task => task
    }
    tasks.size shouldBe 1
    val task = tasks.head
    val standAloneTask = codegen.createStandAloneTask(task)
    val standAloneTaskSource = versionSupport.generateDocument(standAloneTask)

    // update paths of input files - this requires a round-trip de/ser
    // JSON -> IR -> WDL -> update paths -> IR -> JSON
    // which requires some auxilliarly objects
    val outputSerializer: ParameterLinkSerializer = ParameterLinkSerializer(fileResolver, dxApi)
    val taskInputs: Map[DxName, TAT.InputParameter] =
      task.inputs.map(inp => WdlDxName.fromDecodedName(inp.name) -> inp).toMap
    val updatedInputs = inputDeserializer
      .deserializeInputMap(inputs)
      .collect {
        case (dxName, irValue) if !dxName.suffix.exists(_.endsWith(Constants.FlatFilesSuffix)) =>
          val wdlType = taskInputs(dxName).wdlType
          val irType = WdlUtils.toIRType(wdlType)
          val updatedIrValue = Value.transform(irValue, Some(irType), AddBaseDir)
          outputSerializer.createFields(dxName, irType, updatedIrValue)
      }
      .flatten
      .toMap

    // wrap inputs if using manifests
    val finalInputs = updatedInputs match {
      case i if useManifests =>
        Map(
            Constants.InputManifest -> JsObject(i.map {
              case (dxName, jsv) => dxName.decoded -> jsv
            }),
            Constants.OutputId -> JsString("test")
        )
      case i => i
    }

    // create JobMeta
    val jobMeta =
      TaskTestJobMeta(DxWorkerPaths(PosixPath(jobRootDir.toString)),
                      dxApi,
                      logger,
                      finalInputs,
                      instanceTypeDB,
                      standAloneTaskSource,
                      pathToDxFile,
                      useManifests,
                      downloadFiles)

    // create TaskExecutor
    (WdlTaskExecutor.create(jobMeta, streamFiles = streamFiles, checkInstanceType = false), jobMeta)
  }

  // Parse the WDL source code, extract the single task that is supposed to be there,
  // run the task, and compare the outputs to the expected values (if any).
  private def runTask(wdlName: String, useManifests: Boolean = false): TaskTestJobMeta = {
    val (taskExecutor, jobMeta) = createTaskExecutor(wdlName, useManifests = useManifests)
    val outputsExpected = getExpectedOutputs(wdlName)

    // run the steps of task execution in order
    taskExecutor.apply() shouldBe TaskExecutorResult.Success

    if (outputsExpected.isDefined) {
      val outputs = if (useManifests) {
        jobMeta.outputs
          .map(_.get(Constants.OutputManifest) match {
            case Some(value) =>
              val manifestFile = value match {
                case fileObj: JsObject if DxFile.isLinkJson(fileObj) =>
                  DxFile.fromJson(dxApi, fileObj)
                case JsString(uri) if uri.startsWith(DxPath.DxUriPrefix) =>
                  dxApi.resolveFile(uri)
                case other =>
                  throw new Exception(s"invalid manifest file value ${other}")
              }
              // sometimes it takes a while for the output file to close -
              // block here until the file is closed
              if (!Iterator.range(0, 10).exists { i =>
                    if (i > 0) {
                      Thread.sleep(6000)
                    }
                    val desc =
                      dxApi.fileDescribe(manifestFile.id,
                                         Map("fields" -> JsObject("state" -> JsBoolean(true))))
                    desc.fields.get("state") match {
                      case Some(JsString("closed")) => true
                      case _                        => false
                    }
                  }) {
                throw new Exception("manifest file did not close within 60 seconds")
              }
              Manifest
                .parse(
                    new String(jobMeta.dxApi.downloadBytes(manifestFile)).parseJson,
                    WdlDxName
                )
                .jsValues
            case None => Map.empty[DxName, JsValue]
          })
          .getOrElse(Map.empty[DxName, JsValue])
      } else {
        jobMeta.outputs.getOrElse(Map.empty[DxName, JsValue])
      }
      outputs shouldBe outputsExpected.get
    }

    jobMeta
  }

  it should "execute a simple WDL task" in {
    runTask("add")
  }

  it should "execute a WDL task with expressions" in {
    runTask("float_arith")
  }

  it should "evaluate expressions in runtime section" in {
    runTask("expressions_runtime_section")
  }

  it should "evaluate expressions in runtime section II" in {
    runTask("expressions_runtime_section_2")
  }

  it should "evaluate a command section" in {
    runTask("sub")
  }

  it should "run ps in a command section" in {
    runTask("ps")
  }

  it should "localize a file to a task" in {
    runTask("cgrep")
  }

  it should "handle files with same name in different source folders" taggedAs ApiTest in {
    val (taskExecutor, jobMeta) =
      createTaskExecutor("two_files", StreamFiles.None, downloadFiles = false)
    taskExecutor.apply() shouldBe TaskExecutorResult.Success
    // make sure the download manifest has two different folders
    val manifestJs =
      JsUtils.jsFromFile(jobMeta.workerPaths.getDxdaManifestFile().asJavaPath).asJsObject.fields
    manifestJs.size shouldBe 1
    val folders = manifestJs.values.head match {
      case JsArray(array) =>
        array.map {
          case JsObject(file) =>
            file.get("folder") match {
              case Some(JsString(folder)) => folder
              case other                  => throw new Exception(s"invalid manifest entry ${other}")
            }
          case other => throw new Exception(s"invalid manifest entry ${other}")
        }
      case other => throw new Exception(s"invalid manifest ${other}")
    }
    folders.size shouldBe 2
    folders.toSet.size shouldBe 2
  }

  // this test is invalid - automatic coercion to String is not allowed except
  // in string interpolation
  ignore should "handle type coercion" in {
    runTask("cast")
  }

  it should "handle spaces in file paths" in {
    runTask("spaces_in_file_paths")
  }

  it should "read_tsv" in {
    runTask("read_tsv_x")
  }

  it should "write_tsv" in {
    runTask("write_tsv_x")
  }

  it should "optimize task with an empty command section" in {
    val _ = runTask("empty_command_section")
    //task.commandSectionEmpty should be(true)
  }

  it should "handle structs" taggedAs EdgeTest in {
    runTask("Person2")
  }

  it should "handle missing optional files" in {
    runTask("missing_optional_output_file")
  }

  it should "run a python script" in {
    runTask("python_heredoc")
  }

  it should "run a task using manifests" in {
    runTask("add", useManifests = true)
  }

  it should "localize and delocalize files" in {
    val jobMeta = runTask("opt_array")
    jobMeta.outputs shouldBe Some(
        Map(
            WdlDxName.fromSourceName("o") -> JsArray(
                JsObject(
                    "$dnanexus_link" -> JsObject(
                        "id" -> JsString("file-FGqFGBQ0ffPPkYP19gBvFkZy"),
                        "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")
                    )
                ),
                JsObject(
                    "$dnanexus_link" -> JsObject(
                        "id" -> JsString("file-FGqFJ8Q0ffPGVz3zGy4FK02P"),
                        "project" -> JsString("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")
                    )
                )
            )
        )
    )
  }
}
