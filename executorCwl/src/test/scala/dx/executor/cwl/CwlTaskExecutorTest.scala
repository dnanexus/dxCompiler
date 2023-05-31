package dx.executor.cwl

import Assumptions.{cwltoolCallable, dxdaCallable, isLoggedIn}
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
import dx.core.languages.Language
import dx.core.languages.cwl.{CwlDxName, CwlUtils, DxHintSchema}
import dx.cwl.{CommandInputParameter, CommandLineTool, Parser, ParserResult}
import dx.executor.{JobMeta, TaskExecutor, TaskExecutorResult}
import dx.util.protocols.DxFileAccessProtocol
import dx.util.{CodecUtils, FileSourceResolver, FileUtils, JsUtils, Logger, PosixPath, SysUtils}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

import java.nio.file.{Files, Path, Paths}
import scala.util.Random

private case class ToolTestJobMeta(override val workerPaths: DxWorkerPaths,
                                   override val dxApi: DxApi = DxApi.get,
                                   override val logger: Logger = Logger.get,
                                   override val rawJsInputs: Map[DxName, JsValue],
                                   toolName: String,
                                   rawInstanceTypeDb: InstanceTypeDB,
                                   rawSourceCode: String,
                                   pathToDxFile: Map[Path, DxFile] = Map.empty,
                                   useManifestInputs: Boolean = false,
                                   downloadFiles: Boolean = true)
    extends JobMeta(workerPaths, CwlDxName, dxApi, logger) {
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
                                    truncateLogs: Boolean = false,
                                    forwardStd: Boolean = false): Unit = {
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
          // this will never fail due to the way the command script is written - instead
          // we need to read the return code file
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

  override val instanceType: Option[String] = Some(ToolTestJobMeta.InstanceType)

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
      Constants.UseManifests -> JsBoolean(useManifestInputs),
      Constants.PathsAsObjects -> JsTrue,
      Constants.OriginalName -> JsString(toolName)
  )

  override def getExecutableDetail(name: String): Option[JsValue] = {
    executableDetails.get(name)
  }

  override def error(e: Throwable): Unit = {}
}

private object ToolTestJobMeta {
  val InstanceType = "mem_ssd_unicorn"
}

class CwlTaskExecutorTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)
  assume(cwltoolCallable)
  private val logger = Logger.Quiet
  private val dxApi = DxApi()(logger)
  private val unicornInstance = DxInstanceType(
      ToolTestJobMeta.InstanceType,
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
    InstanceTypeDB(Map(ToolTestJobMeta.InstanceType -> unicornInstance))

  private def pathFromBasename(basename: String): Option[Path] = {
    getClass.getResource(s"/tool_executor/${basename}") match {
      case null => None
      case res  => Some(Paths.get(res.getPath))
    }
  }

  object AddBaseDir extends Value.TransformHandler {
    override def apply(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      value match {
        case f: Value.VFile if !f.uri.startsWith(DxPath.DxUriPrefix) =>
          pathFromBasename(f.uri).map(path => f.copy(uri = path.toString)).orElse {
            throw new Exception(s"File ${f.uri} does not exist")
          }
        case f: Value.VFolder if !f.uri.startsWith(DxPath.DxUriPrefix) =>
          pathFromBasename(f.uri).map(path => f.copy(uri = path.toString)).orElse {
            throw new Exception(s"Folder ${f.uri} does not exist")
          }
        case _ => None
      }
    }
  }

  private def getInputs(cwlName: String): Map[DxName, JsValue] = {
    pathFromBasename(s"${cwlName}_input.json") match {
      case Some(path) if Files.exists(path) =>
        JsUtils.getFields(JsUtils.jsFromFile(path)).map {
          case (name, jsv) => CwlDxName.fromEncodedName(name) -> jsv
        }
      case _ => Map.empty
    }
  }

  private def getExpectedOutputs(cwlName: String): Option[Map[DxName, JsValue]] = {
    pathFromBasename(s"${cwlName}_output.json") match {
      case Some(path) if Files.exists(path) =>
        Some(JsUtils.getFields(JsUtils.jsFromFile(path)).map {
          case (name, jsv) => CwlDxName.fromDecodedName(name) -> jsv
        })
      case _ => None
    }
  }

  private def createTaskExecutor(
      cwlName: String,
      useManifests: Boolean,
      pathToDxFile: Map[Path, DxFile]
  ): (CwlTaskExecutor, ToolTestJobMeta) = {
    val cwlFile: Path = pathFromBasename(s"${cwlName}.cwl.json").get
    val inputs = getInputs(cwlName)
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()
    val workerPaths = DxWorkerPaths(PosixPath(jobRootDir.toString))
    workerPaths.createCleanDirs()

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

    val parser = Parser.create(hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClassFromFile(cwlFile) match {
      case (version, Some("CommandLineTool")) if Language.parse(version) == Language.CwlV1_2 => ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL CommandLineTool document of a supported version
               |${cwlFile}""".stripMargin
        )
    }
    val tool = parser.parseFile(cwlFile) match {
      case ParserResult(Some(tool: CommandLineTool), _, _, _) => tool
      case other =>
        throw new Exception(s"expected CWL document to contain a CommandLineTool, not ${other}")
    }

    val taskInputs: Map[DxName, CommandInputParameter] =
      tool.inputs
        .map(inp => CwlDxName.fromSourceName(inp.name) -> inp)
        .toMap
    val outputSerializer: ParameterLinkSerializer =
      ParameterLinkSerializer(fileResolver, dxApi, pathsAsObjects = true)
    val irInputs = inputDeserializer
      .deserializeInputMap(inputs)
      .collect {
        case (dxName, irValue) if !dxName.suffix.exists(_.endsWith(Constants.FlatFilesSuffix)) =>
          val cwlType = taskInputs(dxName).cwlType
          val irType = CwlUtils.toIRType(cwlType)
          val updatedIrValue = Value.transform(irValue, Some(irType), AddBaseDir)
          outputSerializer.createFields(dxName, irType, updatedIrValue)
      }
      .flatten
      .toMap

    // create JobMeta
    val jobMeta =
      ToolTestJobMeta(
          DxWorkerPaths(PosixPath(jobRootDir.toString)),
          dxApi,
          logger,
          irInputs,
          tool.name,
          instanceTypeDB,
          FileUtils.readFileContent(cwlFile),
          pathToDxFile,
          useManifests
      )

    // create TaskExecutor
    (CwlTaskExecutor.create(jobMeta, streamFiles = StreamFiles.None, checkInstanceType = false),
     jobMeta)
  }

  def compareJsv(x: JsValue, y: JsValue, assertEqual: Boolean = true): Int = {
    (x, y) match {
      case (JsObject(fields1), JsObject(fields2)) =>
        val keysCmp = if (assertEqual) {
          withClue("keys") {
            fields1.keys shouldBe fields2.keys
            0
          }
        } else {
          fields1.keys.toVector.sorted
            .zip(fields2.keys.toVector.sorted)
            .iterator
            .map {
              case (k1, k2) => k1.compare(k2)
            }
            .collectFirst {
              case cmp if cmp != 0 => cmp
            }
            .getOrElse(0)
        }
        if (keysCmp != 0) {
          keysCmp
        } else if (fields1.size != fields2.size) {
          fields1.size.compare(fields2.size)
        } else {
          fields1.iterator
            .map {
              case (key, jsv) =>
                withClue(key) {
                  compareJsv(jsv, fields2(key), assertEqual)
                }
            }
            .collectFirst {
              case cmp if cmp != 0 => cmp
            }
            .getOrElse(0)
        }
      case (JsArray(items1), JsArray(items2)) =>
        if (assertEqual) {
          withClue("size") {
            items1.size shouldBe items2.size
          }
        }
        items1
          .sortWith {
            case (j1, j2) => compareJsv(j1, j2, assertEqual = false) < 0
          }
          .zip(items2.sortWith {
            case (j1, j2) => compareJsv(j1, j2, assertEqual = false) < 0
          })
          .iterator
          .zipWithIndex
          .map {
            case ((i1, i2), idx) =>
              withClue(idx) {
                compareJsv(i1, i2, assertEqual)
              }
          }
          .collectFirst {
            case cmp if cmp != 0 => cmp
          }
          .getOrElse(items1.size.compare(items2.size))
      case (JsBoolean(b1), JsBoolean(b2)) if assertEqual =>
        b1 shouldBe b2
        0
      case (JsBoolean(false), JsBoolean(true)) => -1
      case (JsBoolean(true), JsBoolean(false)) => 1
      case (JsBoolean(_), JsBoolean(_))        => 0
      case (JsNumber(n1), JsNumber(n2)) if assertEqual =>
        n1 shouldBe n2
        0
      case (JsNumber(n1), JsNumber(n2)) => n1.compare(n2)
      case (JsString(s1), JsString(s2)) if assertEqual =>
        s1 shouldBe s2
        0
      case (JsString(s1), JsString(s2)) => s1.compare(s2)
      case _ if assertEqual =>
        x shouldBe y
        0
      case _ => x.toString.compare(y.toString)
    }
  }

  private def compareOutputs[T](outputs: Map[T, JsValue], expected: Map[T, JsValue]): Unit = {
    withClue("outputs") {
      outputs.size shouldBe expected.size
      outputs.keys shouldBe expected.keys
      outputs.foreach {
        case (key, value) =>
          withClue(key) {
            compareJsv(value, expected(key))
          }
      }
    }
  }

  // Parse the CWL source code, extract the single tool that is supposed to be there,
  // run the tool, and compare the outputs to the expected values (if any).
  private def runTask(cwlName: String,
                      useManifests: Boolean = false,
                      pathToDxFile: Map[Path, DxFile] = Map.empty): Unit = {
    val (taskExecutor, jobMeta) =
      createTaskExecutor(cwlName, useManifests = useManifests, pathToDxFile = pathToDxFile)
    val outputsExpected = getExpectedOutputs(cwlName)

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
                      Thread.sleep(1000)
                    }
                    val desc =
                      dxApi.fileDescribe(manifestFile.id,
                                         Map("fields" -> JsObject("state" -> JsBoolean(true))))
                    desc.fields.get("state") match {
                      case Some(JsString("closed")) => true
                      case _                        => false
                    }
                  }) {
                throw new Exception("manifest file did not close within 10 seconds")
              }
              Manifest
                .parse(
                    new String(jobMeta.dxApi.downloadBytes(manifestFile)).parseJson,
                    CwlDxName
                )
                .jsValues
            case None => Map.empty[DxName, JsValue]
          })
          .getOrElse(Map.empty[DxName, JsValue])
      } else {
        jobMeta.outputs.getOrElse(Map.empty[DxName, JsValue])
      }
      compareOutputs(outputs, outputsExpected.get)
    }
  }

  it should "handle in-line record type" in {
    runTask("record-in-format")
  }

  it should "handle nested input directory" in {
    runTask(
        "recursive-input-directory",
        pathToDxFile = Map(
            Paths.get("work_dir/a") -> dxApi.file(
                "file-GJ8KZ780yzZyBbf07QVxJY2f",
                Option(DxProject("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")(dxApi))
            ),
            Paths.get("work_dir/b") -> dxApi.file(
                "file-GJ8KZ7j0yzZXkzG57XQkgxbX",
                Option(DxProject("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")(dxApi))
            ),
            Paths.get("work_dir/c/d") -> dxApi.file(
                "file-GJ8KZ800yzZXkzG57XQkgxbZ",
                Option(DxProject("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")(dxApi))
            ),
            Paths.get("work_dir/e") -> dxApi.file(
                "file-GJ8KZ800yzZxbjyX7Q7QgXgV",
                Option(DxProject("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")(dxApi))
            ),
            Paths.get("output.txt") -> dxApi.file(
                "file-GJ713F80yzZy59jJ4p5FY8kG",
                Option(DxProject("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")(dxApi))
            )
        )
    )
  }
}
