package dx.executor

import java.nio.file.{Files, Path}

import dx.{AppException, AppInternalException}
import dx.api.{
  DxAnalysis,
  DxApi,
  DxExecutable,
  DxExecution,
  DxFile,
  DxFileDescCache,
  DxJob,
  DxJobDescribe,
  DxProject,
  DxProjectDescribe,
  Field,
  InstanceTypeDB
}
import dx.core.Constants
import dx.core.io.DxWorkerPaths
import dx.core.ir.Value.VNull
import dx.core.ir.{
  ParameterLink,
  ParameterLinkDeserializer,
  ParameterLinkExec,
  ParameterLinkSerializer,
  Type,
  TypeSerde,
  Value,
  ValueSerde
}
import dx.core.languages.Language
import dx.util.{CodecUtils, FileSourceResolver, FileUtils, JsUtils, Logger, TraceLevel}
import dx.util.protocols.DxFileAccessProtocol
import spray.json._

object JobMeta {
  val inputFile = "job_input.json"
  val outputFile = "job_output.json"
  val errorFile = "job_error.json"
  // this file has all the information about the job that is avaiable
  // at execution time
  val jobInfoFile = "dnanexus-job.json"
  // this file has all the information about the executable
  val executableInfoFile = "dnanexus-executable.json"

  /**
    * Report an error, since this is called from a bash script, we can't simply
    * raise an exception. Instead, we write the error to a standard JSON file.
    * @param e the exception
    */
  def writeError(rootDir: Path, e: Throwable): Unit = {
    val jobErrorPath = rootDir.resolve(errorFile)
    val errType = e match {
      case _: AppException         => "AppError"
      case _: AppInternalException => "AppInternalError"
      case _: Throwable            => "AppInternalError"
    }
    // We are limited in what characters can be written to json, so we
    // provide a short description for json.
    //
    // Note: we sanitize this string, to be absolutely sure that
    // it does not contain problematic JSON characters.
    val errMsg = JsObject(
        "error" -> JsObject(
            "type" -> JsString(errType),
            "message" -> JsUtils.sanitizedString(e.getMessage)
        )
    ).prettyPrint
    FileUtils.writeFileContent(jobErrorPath, errMsg)
  }
}

/**
  * Encapsulates all the metadata used by the executor, including:
  * 1. metadata files on the worker
  * 2. job details
  * 3. application details
  * @param workerPaths DxWorkerPaths
  * @param dxApi DxApi
  * @param logger Logger
  */
abstract class JobMeta(val workerPaths: DxWorkerPaths, val dxApi: DxApi, val logger: Logger) {
  def project: DxProject

  lazy val projectDesc: DxProjectDescribe = project.describe()

  def jsInputs: Map[String, JsValue]

  private lazy val dxFileDescCache: DxFileDescCache = {
    val allFilesReferenced = jsInputs.flatMap {
      case (_, jsElem) => DxFile.findFiles(dxApi, jsElem)
    }.toVector
    // Describe all the files and build a lookup cache
    DxFileDescCache(dxApi.describeFilesBulk(allFilesReferenced))
  }

  lazy val inputDeserializer: ParameterLinkDeserializer =
    ParameterLinkDeserializer(dxFileDescCache, dxApi)
  lazy val inputs: Map[String, Value] = {
    // if we have access to the inputSpec, use it to guide deserialization
    getExecutableAttribute("inputSpec") match {
      case Some(JsArray(spec)) =>
        val inputTypes = TypeSerde.fromNativeSpec(spec)
        jsInputs.map {
          case (key, value) =>
            inputTypes.get(key) match {
              case None =>
                logger.warning(s"inputSpec is missing field ${key}")
                key -> inputDeserializer.deserializeInput(value)
              case Some(Type.THash) =>
                logger.trace(
                    s"""expected type of input field ${key} is THash, which may represent an
                       |unknown schema type, so deserializing without type""".stripMargin
                      .replaceAll("\n", " ")
                )
                key -> inputDeserializer.deserializeInput(value)
              case Some(t) =>
                key -> inputDeserializer.deserializeInputWithType(value, t)
            }
        }
      case None =>
        inputDeserializer.deserializeInputMap(jsInputs)
      case other =>
        throw new Exception(s"invalid inputSpec ${other}")
    }
  }

  lazy val fileResolver: FileSourceResolver = {
    val dxProtocol = DxFileAccessProtocol(dxApi, dxFileDescCache)
    val fileResolver = FileSourceResolver.create(
        localDirectories = Vector(workerPaths.getWorkDir()),
        userProtocols = Vector(dxProtocol),
        logger = logger
    )
    FileSourceResolver.set(fileResolver)
    fileResolver
  }

  def writeJsOutputs(outputJs: Map[String, JsValue]): Unit

  lazy val outputSerializer: ParameterLinkSerializer = ParameterLinkSerializer(fileResolver, dxApi)

  private def outputTypesEqual(expectedType: Type, actualType: Type): Boolean = {
    val comparisonType = if (Type.isNative(actualType)) {
      actualType
    } else {
      Type.THash
    }
    comparisonType == expectedType || (
        Type.isOptional(expectedType) && Type.ensureOptional(comparisonType) == expectedType
    )
  }

  def writeOutputs(outputs: Map[String, (Type, Value)]): Unit = {
    getExecutableAttribute("outputSpec").foreach {
      case JsArray(spec) =>
        // check that the actual types are the same as the expected types
        val outputTypes = TypeSerde.fromNativeSpec(spec)
        outputs.foreach {
          case (key, (actualType, _)) =>
            outputTypes.get(key) match {
              case Some(t) if outputTypesEqual(t, actualType) => ()
              case Some(t) if t == Type.THash =>
                logger.trace(
                    s"""expected type of output field ${key} is THash which may represent an
                       |unknown schema type, so deserializing without type""".stripMargin
                      .replaceAll("\n", " ")
                )
              case Some(t) =>
                throw new Exception(
                    s"""output field ${key} has mismatch between actual type ${actualType}
                       |and expected type ${t}""".stripMargin.replaceAll("\n", " ")
                )
              case None =>
                logger.warning(s"outputSpec is missing field ${key}")
            }
        }
      case other =>
        throw new Exception(s"invalid outputSpec ${other}")
    }
    // write outputs, ignore null values - these could occur for optional
    // values that were not specified.
    val outputJs = outputs
      .collect {
        case (name, (t, v)) if v != VNull =>
          outputSerializer.createFields(name, t, v)
      }
      .flatten
      .toMap
    writeJsOutputs(outputJs)
  }

  def createOutputLinks(outputs: Map[String, (Type, Value)]): Map[String, ParameterLink] = {
    outputs.collect {
      case (name, (t, v)) if v != VNull =>
        name -> outputSerializer.createLink(t, v)
    }
  }

  def serializeOutputLinks(outputs: Map[String, ParameterLink]): Map[String, JsValue] = {
    outputs
      .map {
        case (name, link) => outputSerializer.createFieldsFromLink(link, name)
      }
      .flatten
      .toMap
  }

  def writeOutputLinks(outputs: Map[String, ParameterLink]): Unit = {
    writeJsOutputs(serializeOutputLinks(outputs))
  }

  def createOutputLinks(execution: DxExecution,
                        irOutputFields: Map[String, Type],
                        prefix: Option[String] = None): Map[String, ParameterLink] = {
    irOutputFields.map {
      case (name, t) =>
        val fqn = prefix.map(p => s"${p}.${name}").getOrElse(name)
        fqn -> ParameterLinkExec(execution, name, t)
    }
  }

  def serializeOutputLinks(execution: DxExecution,
                           irOutputFields: Map[String, Type]): Map[String, JsValue] = {
    createOutputLinks(execution, irOutputFields)
      .flatMap {
        case (name, link) => outputSerializer.createFieldsFromLink(link, name)
      }
      .filter {
        case (_, null | JsNull) => false
        case _                  => true
      }
  }

  def writeOutputLinks(execution: DxExecution, irOutputFields: Map[String, Type]): Unit = {
    writeJsOutputs(serializeOutputLinks(execution, irOutputFields))
  }

  def jobId: String

  def analysis: Option[DxAnalysis]

  def parentJob: Option[DxJob]

  def instanceType: Option[String]

  def getJobDetail(name: String): Option[JsValue]

  def getExecutableAttribute(name: String): Option[JsValue]

  def getExecutableDetail(name: String): Option[JsValue]

  lazy val language: Option[Language.Language] = getExecutableDetail(Constants.Language) match {
    case Some(JsString(lang)) => Some(Language.withName(lang))
    case None =>
      logger.warning("This applet ws built with an old version of dxWDL - please rebuild")
      // we will attempt to detect the language/version later
      None
    case other =>
      throw new Exception(s"unexpected ${Constants.Language} value ${other}")
  }

  lazy val sourceCode: String = {
    val sourceCodeEncoded = getExecutableDetail(Constants.SourceCode) match {
      case Some(JsString(s)) => s
      case None =>
        logger.warning("This applet ws built with an old version of dxWDL - please rebuild")
        val JsString(s) =
          getExecutableDetail("wdlSourceCode")
            .orElse(getExecutableDetail("womSourceCode"))
            .getOrElse(
                throw new Exception("executable details does not contain source code")
            )
        s
      case other =>
        throw new Exception(s"unexpected ${Constants.SourceCode} value ${other}")
    }
    CodecUtils.base64DecodeAndGunzip(sourceCodeEncoded)
  }

  lazy val instanceTypeDb: InstanceTypeDB = getExecutableDetail(Constants.InstanceTypeDb) match {
    case Some(JsString(s)) =>
      val js = CodecUtils.base64DecodeAndGunzip(s)
      js.parseJson.convertTo[InstanceTypeDB]
    case other =>
      throw new Exception(s"unexpected ${Constants.InstanceTypeDb} value ${other}")
  }

  lazy val defaultRuntimeAttrs: Map[String, Value] =
    getExecutableDetail(Constants.RuntimeAttributes) match {
      case Some(JsObject(fields)) => ValueSerde.deserializeMap(fields)
      case Some(JsNull) | None    => Map.empty
      case other =>
        throw new Exception(s"unexpected ${Constants.RuntimeAttributes} value ${other}")
    }

  lazy val delayWorkspaceDestruction: Option[Boolean] =
    getExecutableDetail("delayWorkspaceDestruction") match {
      case Some(JsBoolean(flag)) => Some(flag)
      case None                  => None
    }

  lazy val blockPath: Vector[Int] = getExecutableDetail("blockPath") match {
    case Some(JsArray(arr)) if arr.nonEmpty =>
      arr.map {
        case JsNumber(n) => n.toInt
        case _           => throw new Exception("Bad value ${arr}")
      }
    case None =>
      Vector.empty
    case other =>
      throw new Exception(s"Bad value ${other}")
  }

  lazy val scatterStart: Int = getJobDetail(Constants.ContinueStart) match {
    case Some(JsNumber(s)) => s.toIntExact
    case Some(other) =>
      throw new Exception(s"Invalid value ${other} for  ${Constants.ContinueStart}")
    case _ => 0
  }

  lazy val scatterSize: Int = getExecutableDetail(Constants.ScatterChunkSize) match {
    case Some(JsNumber(n)) => n.toIntExact
    case None              => Constants.JobPerScatterDefault
    case other =>
      throw new Exception(s"Invalid value ${other} for ${Constants.ScatterChunkSize}")
  }

  def error(e: Throwable): Unit
}

case class WorkerJobMeta(override val workerPaths: DxWorkerPaths = DxWorkerPaths.default,
                         override val dxApi: DxApi = DxApi.get,
                         override val logger: Logger = Logger.get)
    extends JobMeta(workerPaths, dxApi, logger) {
  lazy val project: DxProject = dxApi.currentProject

  private val rootDir = workerPaths.getRootDir()
  private val inputPath = rootDir.resolve(JobMeta.inputFile)
  private val outputPath = rootDir.resolve(JobMeta.outputFile)
  private val jobInfoPath = rootDir.resolve(JobMeta.jobInfoFile)
  private val executableInfoPath = rootDir.resolve(JobMeta.executableInfoFile)

  lazy val jsInputs: Map[String, JsValue] = {
    if (Files.exists(inputPath)) {
      JsUtils.getFields(JsUtils.jsFromFile(inputPath))
    } else {
      logger.warning(s"input meta-file ${inputPath} does not exist")
      Map.empty
    }
  }

  def writeJsOutputs(outputJs: Map[String, JsValue]): Unit = {
    JsUtils.jsToFile(JsObject(outputJs), outputPath)
  }

  private lazy val jobInfo: Map[String, JsValue] = {
    if (Files.exists(jobInfoPath)) {
      FileUtils.readFileContent(jobInfoPath).parseJson.asJsObject.fields
    } else {
      logger.warning(s"info meta-file ${jobInfoPath} does not exist")
      Map.empty
    }
  }

  def jobId: String = DxApi.get.currentJob.id

  private lazy val jobDesc: DxJobDescribe = DxApi.get.currentJob.describe(
      Set(Field.Executable, Field.Details, Field.InstanceType)
  )

  private lazy val jobDetails: Map[String, JsValue] =
    jobDesc.details.map(_.asJsObject.fields).getOrElse(Map.empty)

  private lazy val executable: DxExecutable = jobInfo.get("executable") match {
    case None =>
      logger.trace(
          "executable field not found locally, performing an API call.",
          minLevel = TraceLevel.None
      )
      jobDesc.executable.get
    case Some(JsString(x)) if x.startsWith("app-") =>
      dxApi.app(x)
    case Some(JsString(x)) if x.startsWith("applet-") =>
      dxApi.applet(x)
    case Some(other) =>
      throw new Exception(s"Malformed executable field ${other} in job info")
  }

  private lazy val executableInfo: Map[String, JsValue] = {
    if (Files.exists(executableInfoPath)) {
      FileUtils.readFileContent(executableInfoPath).parseJson.asJsObject.fields
    } else {
      logger.warning(s"executable meta-file ${executableInfoPath} does not exist")
      Map.empty
    }
  }

  override def getExecutableAttribute(name: String): Option[JsValue] = {
    executableInfo.get(name)
  }

  private lazy val executableDetails: Map[String, JsValue] = {
    executable.describe(Set(Field.Details)).details.get.asJsObject.fields
  }

  def analysis: Option[DxAnalysis] = jobDesc.analysis

  def parentJob: Option[DxJob] = jobDesc.parentJob

  def instanceType: Option[String] = jobDesc.instanceType

  def getJobDetail(name: String): Option[JsValue] = jobDetails.get(name)

  def getExecutableDetail(name: String): Option[JsValue] = executableDetails.get(name)

  def error(e: Throwable): Unit = {
    JobMeta.writeError(rootDir, e)
  }
}
