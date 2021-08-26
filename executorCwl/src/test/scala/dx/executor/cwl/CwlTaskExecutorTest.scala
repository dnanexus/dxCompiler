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
import dx.cwl.{CommandInputParameter, CommandLineTool, Parser}
import dx.executor.{JobMeta, TaskExecutor}
import dx.util.protocols.DxFileAccessProtocol
import dx.util.{CodecUtils, FileSourceResolver, FileUtils, JsUtils, Logger, SysUtils}
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
                                   useManifestInputs: Boolean = false)
    extends JobMeta(workerPaths, dxApi, logger) {
  var outputs: Option[Map[DxName, JsValue]] = None

  override val project: DxProject = null

  override def writeRawJsOutputs(outputJs: Map[DxName, JsValue]): Unit = {
    outputs = Some(outputJs)
  }

  override lazy val jobId: String = s"job-${Random.alphanumeric.take(24).mkString}"

  override def runJobScriptFunction(name: String,
                                    successCodes: Option[Set[Int]] = Some(Set(0))): Unit = {
    name match {
      case TaskExecutor.DownloadDxda if dxdaCallable =>
        val dxdaManifest = workerPaths.getDxdaManifestFile().toString
        SysUtils.execCommand(
            s"""bzip2 ${dxdaManifest} &&
               |cd / &&
               |dx-download-agent download ${dxdaManifest}.bz2""".stripMargin
              .replaceAll("\n", " ")
        )
      case TaskExecutor.RunCommand =>
        val script: Path = workerPaths.getCommandFile()
        if (Files.exists(script)) {
          // this will never fail due to the way the command script is written - instead
          // we need to read the return code file
          val (_, stdout, stderr) =
            SysUtils.execCommand(script.toString, exceptionOnFailure = false)
          val rc = FileUtils.readFileContent(workerPaths.getReturnCodeFile()).trim.toInt
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

  override def getExecutableAttribute(name: String): Option[JsValue] = {
    name match {
      case "name" => Some(JsString(toolName))
      case _      => None
    }
  }

  private val executableDetails: Map[String, JsValue] = Map(
      Constants.InstanceTypeDb -> JsString(
          CodecUtils.gzipAndBase64Encode(
              rawInstanceTypeDb.toJson.prettyPrint
          )
      ),
      Constants.SourceCode -> JsString(CodecUtils.gzipAndBase64Encode(rawSourceCode)),
      Constants.UseManifests -> JsBoolean(useManifestInputs),
      Constants.PathsAsObjects -> JsTrue
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
      Some(1.00f)
  )
  private val instanceTypeDB =
    InstanceTypeDB(Map(ToolTestJobMeta.InstanceType -> unicornInstance), pricingAvailable = true)

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
          case (name, jsv) => CwlDxName.fromEncodedParameterName(name) -> jsv
        }
      case _ => Map.empty
    }
  }

  private def getExpectedOutputs(cwlName: String): Option[Map[String, JsValue]] = {
    pathFromBasename(s"${cwlName}_output.json") match {
      case Some(path) if Files.exists(path) => Some(JsUtils.getFields(JsUtils.jsFromFile(path)))
      case _                                => None
    }
  }

  private def createTaskExecutor(
      cwlName: String,
      useManifests: Boolean,
      streamFiles: StreamFiles.StreamFiles = StreamFiles.None,
      waitOnUpload: Boolean = true
  ): (CwlTaskExecutor, ToolTestJobMeta) = {
    val cwlFile: Path = pathFromBasename(s"${cwlName}.cwl").get
    val inputs = getInputs(cwlName)
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()
    val workerPaths = DxWorkerPaths(jobRootDir)
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
    parser.detectVersionAndClass(cwlFile) match {
      case Some((version, "CommandLineTool")) if Language.parse(version) == Language.CwlV1_2 => ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL CommandLineTool document of a supported version
               |${cwlFile}""".stripMargin
        )
    }
    val tool = parser.parseFile(cwlFile) match {
      case (tool: CommandLineTool, _) => tool
      case other =>
        throw new Exception(s"expected CWL document to contain a CommandLineTool, not ${other}")
    }

    val taskInputs: Map[DxName, CommandInputParameter] =
      tool.inputs
        .map(inp => CwlDxName.fromRawParameterName(inp.name) -> inp)
        .toMap
    val outputSerializer: ParameterLinkSerializer =
      ParameterLinkSerializer(fileResolver, dxApi, pathsAsObjects = true)
    val irInputs = inputDeserializer
      .deserializeInputMap(inputs)
      .collect {
        case (dxName, irValue) if !dxName.suffix.contains(Constants.FlatFilesSuffix) =>
          val cwlType = taskInputs(dxName).cwlType
          val irType = CwlUtils.toIRType(cwlType)
          val updatedIrValue = Value.transform(irValue, Some(irType), AddBaseDir)
          outputSerializer.createFields(dxName, irType, updatedIrValue)
      }
      .flatten
      .toMap

    // create JobMeta
    val jobMeta =
      ToolTestJobMeta(DxWorkerPaths(jobRootDir),
                      dxApi,
                      logger,
                      irInputs,
                      tool.name,
                      instanceTypeDB,
                      FileUtils.readFileContent(cwlFile),
                      useManifests)

    // create TaskExecutor
    (CwlTaskExecutor.create(jobMeta,
                            streamFiles = streamFiles,
                            waitOnUpload = waitOnUpload,
                            checkInstanceType = false),
     jobMeta)
  }

  // Parse the CWL source code, extract the single tool that is supposed to be there,
  // run the tool, and compare the outputs to the expected values (if any).
  private def runTask(cwlName: String, useManifests: Boolean = false): Unit = {
    val (taskExecutor, jobMeta) = createTaskExecutor(cwlName, useManifests = useManifests)
    val outputsExpected = getExpectedOutputs(cwlName)

    // run the steps of task execution in order
    taskExecutor.apply() shouldBe true

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
              val manifest =
                Manifest.parse(new String(jobMeta.dxApi.downloadBytes(manifestFile)).parseJson)
              manifest.jsValues
            case None => Map.empty[String, JsValue]
          })
          .getOrElse(Map.empty[String, JsValue])
      } else {
        jobMeta.outputs.getOrElse(Map.empty[String, JsValue])
      }
      outputs shouldBe outputsExpected.get
    }
  }

  it should "handle in-line record type" in {
    runTask("record-in-format")
  }

  it should "handle nested input directory" in {
    runTask("recursive-input-directory")
  }
}
