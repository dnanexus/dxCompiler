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
  DxPath,
  DxProject,
  DxProjectDescribe,
  Field,
  InstanceTypeDB
}
import dx.core.Constants
import dx.core.io.DxWorkerPaths
import dx.core.ir.Value.VNull
import dx.core.ir.{
  ExecutableLink,
  Manifest,
  Parameter,
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

import scala.annotation.tailrec

object JobMeta {
  val InputFile = "job_input.json"
  val OutputFile = "job_output.json"
  val ErrorFile = "job_error.json"
  // this file has all the information about the job that is avaiable
  // at execution time
  val JobInfoFile = "dnanexus-job.json"
  // this file has all the information about the executable
  val ExecutableInfoFile = "dnanexus-executable.json"

  /**
    * Report an error, since this is called from a bash script, we can't simply
    * raise an exception. Instead, we write the error to a standard JSON file.
    * @param rootDir the directory that contains the error file
    * @param e the exception
    */
  def writeError(rootDir: Path, e: Throwable): Unit = {
    val jobErrorPath = rootDir.resolve(ErrorFile)
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

  def jobId: String

  lazy val manifestFolder = s"${dxApi.currentProject.id}:/.d/${jobId}"

  def rawJsInputs: Map[String, JsValue]

  /**
    * If the task/workflow was compiled with -useManifests, then the job inputs
    * will be manifest files, rather than the actual inputs expected by the
    * task/workflow. The inputs will be:
    * 1. At least one manifest file and/or a manifest as a JSON hash
    * 2. A linking hash (required if there is more than one manifest)
    * 3. Zero or more workflow manifest files and/or a workflow manifest as a
    *    JSON hash
    * 4. A linking hash for the workflow manifest files (required if there is
    *    more than one workflow manifest)
    * We need to:
    * 1. Download and parse all the manifest files
    * 2. If there is no linking hash, then just return the values from the
    *    single manifest
    * 3. Otherwise, resolve all the links. A link maps a task/workflow input
    *    name to a value. The link value is one of the following:
    *    * { "value\_\_\_": <value> } - a constant value
    *    * { "<manifest_id>": "<manifest_param_name>" } - a reference to a
    *      parameter that is looked up in the manifest(s)
    *    * { "workflow\_\_\_": "<workflow_param_name>" } - a reference to a
    *      workflow input parameter that is looked up in the workflow
    *      manifest(s)
    * All missing parameters are ignored since they may be optional.
    * @param rawInputs raw job inputs
    * @return the actual task/workflow inputs extracted from the manifest
    */
  private def unpackManifests(rawInputs: Map[String, JsValue]): Map[String, JsValue] = {
    logger.trace(s"""Unpacking manifests from raw inputs:
                    |${rawInputs}""".stripMargin)

    // resolve a file input to a DxFile
    def resolveManifestFiles(key: String): Vector[DxFile] = {
      rawInputs.get(key) match {
        case Some(JsArray(files)) =>
          files.map {
            case fileObj: JsObject if DxFile.isLinkJson(fileObj) =>
              DxFile.fromJson(dxApi, fileObj)
            case JsString(uri) if uri.startsWith(DxPath.DxUriPrefix) =>
              dxApi.resolveFile(uri)
            case other =>
              throw new Exception(s"invalid manifest file value ${other}")
          }
        case None => Vector.empty
        case other =>
          throw new Exception(s"invalid value ${other} for ${key} field")
      }
    }

    val manifestHash = rawInputs.get(Constants.InputManifest) match {
      case Some(hash: JsObject) => Some(hash)
      case None                 => None
      case other =>
        throw new Exception(s"invalid manifest ${other}")
    }
    val manifestFiles = resolveManifestFiles(Constants.InputManifestFiles)
    // links have been serialized, so they actually look like:
    // { "___": { "x": { "wrapped___": { "value___": ... } } }
    // we do not want to deserialize the values (because we don't have the
    // type information here), so instead we have to manually unwrap the fields
    val manifestLinks = rawInputs.get(Constants.InputLinks) match {
      case Some(JsObject(fields)) if fields.contains(Parameter.ComplexValueKey) =>
        fields(Parameter.ComplexValueKey).asJsObject.fields
      case Some(JsObject(fields)) => fields
      case None                   => Map.empty
      case other =>
        throw new Exception(s"invalid manifest links ${other}")
    }

    if (manifestLinks.isEmpty) {
      // if there are no manifest links, then there must be at most one manifest
      // (either as a hash or a file). We can just parse it and return the values.
      val manifestJson = (manifestHash, manifestFiles) match {
        case (None, Vector())         => None
        case (Some(inputs), Vector()) => Some(inputs)
        case (None, Vector(manifestFile)) =>
          Some(new String(dxApi.downloadBytes(manifestFile)).parseJson)
        case _ =>
          throw new Exception(
              "manifest links are required when there is more than one manifest file"
          )
      }
      return manifestJson.map(Manifest.parse(_).jsValues).getOrElse(Map.empty)
    }

    def downloadAndParseManifests(files: Vector[DxFile]): Vector[Manifest] = {
      // download manifest bytes, convert to JSON, and parse into Manifest object
      val manifests = files
        .map(dxFile => Manifest.parse(new String(dxApi.downloadBytes(dxFile)).parseJson))
      // if there is more than one manifest, check that they all have IDs
      if (manifests.size > 1) {
        manifests.foreach {
          case manifest: Manifest if manifest.id.isDefined => ()
          case _ =>
            throw new Exception("when there are multiple manifests, all manifests must have ID")
        }
      }
      manifests
    }

    // parse all the manifests and combine into a single Vector
    val manifests = manifestHash.map(h => Vector(Manifest.parse(h))).getOrElse(Vector.empty) ++
      downloadAndParseManifests(manifestFiles)
    // create lookup function for manifests
    val manifestLookup: (String, String) => Option[JsValue] = if (manifests.size == 1) {
      val values = manifests.head.jsValues
      (_: String, manifestParamName: String) => values.get(manifestParamName)
    } else {
      val manifestMap = manifests.map {
        case m: Manifest if m.id.isDefined => m.id.get -> m.jsValues
        case _ =>
          throw new Exception(
              "when there are multiple manifests, all manifests must have id defined"
          )
      }.toMap
      (manifestId: String, manifestParamName: String) => {
        manifestMap.get(manifestId).flatMap(_.get(manifestParamName))
      }
    }

    // create lookup function for workflow manifests - this is lazy because the
    // manifest linkes may not contain any workflow references
    lazy val workflowManifestLookup: String => Option[JsValue] = {
      val workflowManifestHash = rawInputs.get(Constants.WorkflowInputManifest) match {
        case Some(hash: JsObject) => Some(hash)
        case None                 => None
        case other =>
          throw new Exception(s"invalid manifest ${other}")
      }
      val workflowManifestFiles = resolveManifestFiles(Constants.WorkflowInputManifestFiles)
      val workflowManifests =
        workflowManifestHash.map(h => Vector(Manifest.parse(h))).getOrElse(Vector.empty) ++
          downloadAndParseManifests(workflowManifestFiles)
      val workflowManifestLinks = rawInputs.get(Constants.WorkflowInputLinks) match {
        case Some(JsObject(fields)) if fields.contains(Parameter.ComplexValueKey) =>
          fields(Parameter.ComplexValueKey).asJsObject.fields
        case Some(JsObject(fields)) => fields
        case None                   => Map.empty[String, JsValue]
        case other =>
          throw new Exception(s"invalid workflow manifest links ${other}")
      }
      if (workflowManifests.isEmpty) {
        throw new Exception("there are no workflow manifest files")
      } else if (workflowManifestLinks.isEmpty) {
        if (workflowManifests.size == 1) {
          val values = workflowManifests.head.jsValues
          (paramName: String) => values.get(paramName)
        } else {
          throw new Exception(
              "when there are multiple workflow manifests, workflow manifest links are required"
          )
        }
      } else {
        val workflowManifestMap = workflowManifests.map {
          case m: Manifest if m.id.isDefined => m.id.get -> m.jsValues
          case _ =>
            throw new Exception(
                "when there are multiple workflow manifests, all manifests must have id defined"
            )
        }.toMap
        @tailrec
        def workflowManifestLookup(paramName: String): Option[JsValue] = {
          workflowManifestLinks.get(paramName) match {
            case Some(JsObject(fields)) if fields.keySet == Set(ValueSerde.WrappedValueKey) =>
              fields(ValueSerde.WrappedValueKey).asJsObject.fields.head match {
                case (Constants.ValueKey, value) => Some(value)
                case (Constants.WorkflowKey, JsString(workflowParamName)) =>
                  workflowManifestLookup(workflowParamName)
                case (manifestId, JsString(manifestParamName)) =>
                  workflowManifestMap.get(manifestId).flatMap(_.get(manifestParamName))
                case other =>
                  throw new Exception(s"invalid workflow manifest link entry ${other}")
              }
            case other =>
              throw new Exception(s"invalid manifest link ${other}")
          }
        }
        workflowManifestLookup
      }
    }

    // lookup all values in the manifest links
    manifestLinks.flatMap {
      case (paramName, JsObject(fields)) if fields.keySet == Set(ValueSerde.WrappedValueKey) =>
        val value = fields(ValueSerde.WrappedValueKey).asJsObject.fields.head match {
          case (Constants.ValueKey, value) => Some(value)
          case (Constants.WorkflowKey, JsString(workflowParamName)) =>
            workflowManifestLookup(workflowParamName)
          case (manifestId, JsString(manifestParamName)) =>
            manifestLookup(manifestId, manifestParamName)
          case other =>
            throw new Exception(s"invalid manifest link entry ${other}")
        }
        value.map(paramName -> _)
      case other =>
        throw new Exception(s"invalid manifest link ${other}")
    }
  }

  lazy val jsInputs: Map[String, JsValue] = {
    if (useManifests) {
      unpackManifests(rawJsInputs)
    } else {
      rawJsInputs
    }
  }

  protected lazy val dxFileDescCache: DxFileDescCache = {
    val allFilesReferenced = jsInputs.flatMap {
      case (_, jsElem) => DxFile.findFiles(dxApi, jsElem)
    }.toVector
    // Describe all the files and build a lookup cache
    DxFileDescCache(dxApi.describeFilesBulk(allFilesReferenced))
  }

  lazy val inputDeserializer: ParameterLinkDeserializer =
    ParameterLinkDeserializer(dxFileDescCache, dxApi)

  lazy val inputSpec: Map[String, Type] = {
    getExecutableAttribute("inputSpec")
      .map {
        case JsArray(spec) => TypeSerde.fromNativeSpec(spec)
        case other         => throw new Exception(s"invalid inputSpec ${other}")
      }
      .getOrElse(Map.empty)
  }

  lazy val inputs: Map[String, Value] = {
    // If we are using manifests, then the inputSpec won't match the task inputs.
    // Otherwise, if we have access to the inputSpec, use it to guide deserialization.
    jsInputs.map {
      case (key, value) if useManifests =>
        key -> inputDeserializer.deserializeInput(value)
      case (key, value) =>
        val irValue = inputSpec.get(key) match {
          case None =>
            logger.warning(s"inputSpec is missing field ${key}")
            inputDeserializer.deserializeInput(value)
          case Some(Type.THash) =>
            logger.trace(
                s"""expected type of input field ${key} is THash, which may represent an
                   |unknown schema type, so deserializing without type""".stripMargin
                  .replaceAll("\n", " ")
            )
            inputDeserializer.deserializeInput(value)
          case Some(t) =>
            inputDeserializer.deserializeInputWithType(value, t)
        }
        key -> irValue
    }
  }

  lazy val primaryInputs: Map[String, Value] = {
    // discard auxiliary fields
    inputs.view.filterKeys(!_.endsWith(ParameterLink.FlatFilesSuffix)).toMap
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

  protected def writeRawJsOutputs(outputJs: Map[String, JsValue]): Unit

  def writeJsOutputs(outputJs: Map[String, JsValue]): Unit = {
    val rawOutputJs = if (useManifests && !outputJs.contains(Constants.OutputManifest)) {
      val manifestId = rawJsInputs.get(Constants.OutputId) match {
        case Some(JsString(id)) => id
        case other =>
          throw new Exception(s"missing or invalid outputId ${other}")
      }
      val manifestValues = rawJsInputs.get(Constants.CallName) match {
        case Some(JsString(callName)) =>
          val prefix = s"${callName}${Parameter.ComplexValueKey}"
          outputJs.map {
            case (name, value) => s"${prefix}${name}" -> value
          }
        case None  => outputJs
        case other => throw new Exception(s"invalid call name value ${other}")
      }
      val manifest = Manifest(manifestValues, id = Some(manifestId))
      val destination = s"${manifestFolder}/${jobId}_output.manifest.json"
      val manifestDxFile = dxApi.uploadString(manifest.toJson.prettyPrint, destination)
      outputSerializer
        .createFields(Constants.OutputManifest, Type.TFile, Value.VFile(manifestDxFile.asUri))
        .toMap
    } else {
      outputJs
    }
    writeRawJsOutputs(rawOutputJs)
  }

  lazy val outputSerializer: ParameterLinkSerializer = ParameterLinkSerializer(fileResolver, dxApi)

  /**
    * Prepares the inputs for a subjob. If using manifests, creates a manifest
    * with the same output ID as the current job.
    * @param inputs raw subjob inputs
    * @param executableLink the subjob executable
    * @param callName the call name to use to prefix the subjob outputs
    * @return serialized subjob inputs
    */
  def prepareSubjobInputs(inputs: Map[String, (Type, Value)],
                          executableLink: ExecutableLink,
                          callName: Option[String] = None): Map[String, JsValue] = {
    val inputsJs = outputSerializer.createFieldsFromMap(inputs)
    // Check that we have all the compulsory arguments.
    // Note that we don't have the information here to tell difference between optional and non-
    // optional parameters, so we emit warnings for missing arguments.
    executableLink.inputs.keys.filterNot(inputsJs.contains).foreach { argName =>
      logger.warning(s"Missing argument ${argName} to call ${executableLink.name}", force = true)
    }
    if (useManifests) {
      val requiredInputs = Map(
          Constants.InputManifest -> JsObject(inputsJs),
          Constants.OutputId -> rawJsInputs(Constants.OutputId)
      )
      val callNameInputs =
        callName.map(name => Map(Constants.CallName -> JsString(name))).getOrElse(Map.empty)
      requiredInputs ++ callNameInputs
    } else {
      inputsJs
    }
  }

  @tailrec
  private def outputTypesEqual(expectedType: Type, actualType: Type): Boolean = {
    (expectedType, actualType) match {
      case (a, b) if a == b       => true
      case (a: Type.TOptional, b) =>
        // non-optional actual type is compatible with optional expected type
        outputTypesEqual(Type.unwrapOptional(a), Type.unwrapOptional(b))
      case (Type.TArray(a, false), Type.TArray(b, true)) =>
        // non-empty actual array is compatible with maybe-empty expected array
        outputTypesEqual(a, b)
      case (a, b) if !Type.isNative(b) =>
        // non-native types are represented as hashes
        outputTypesEqual(a, Type.THash)
      case _ => false
    }
  }

  lazy val outputSpec: Map[String, Type] = {
    getExecutableAttribute("outputSpec")
      .map {
        case JsArray(spec) => TypeSerde.fromNativeSpec(spec)
        case other         => throw new Exception(s"invalid outputSpec ${other}")
      }
      .getOrElse(Map.empty)
  }

  private def validateOutput(name: String, actualType: Type): Unit = {
    outputSpec.get(name) match {
      case Some(t) if outputTypesEqual(t, actualType) => ()
      case Some(t) if t == Type.THash =>
        logger.trace(
            s"""expected type of output field ${name} is THash which may represent an
               |unknown schema type, so deserializing without type""".stripMargin
              .replaceAll("\n", " ")
        )
      case Some(t) =>
        throw new Exception(
            s"""output field ${name} has mismatch between actual type ${actualType}
               |and expected type ${t}""".stripMargin.replaceAll("\n", " ")
        )
      case None =>
        logger.warning(s"outputSpec is missing field ${name}")
    }
  }

  def writeOutputs(outputs: Map[String, (Type, Value)]): Unit = {
    // If we are using manifests, then the outputSpec won't match the task outputs.
    // Otherwise, if we have access to the outputSpec, use it to validate the outputs.
    if (!useManifests) {
      outputs.foreach {
        case (name, (actualType, _)) => validateOutput(name, actualType)
      }
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

  def createOutputLinks(outputs: Map[String, (Type, Value)],
                        validate: Boolean = true): Map[String, ParameterLink] = {
    outputs.collect {
      case (name, (actualType, value)) if value != VNull =>
        if (validate) {
          validateOutput(name, actualType)
        }
        name -> outputSerializer.createLink(actualType, value)
    }
  }

  def writeOutputLinks(outputs: Map[String, ParameterLink]): Unit = {
    val serializedLinks = outputs
      .map {
        case (name, link) => outputSerializer.createFieldsFromLink(link, name)
      }
      .flatten
      .toMap
    writeJsOutputs(serializedLinks)
  }

  def createExecutionOutputLinks(execution: DxExecution,
                                 irOutputFields: Map[String, Type],
                                 prefix: Option[String] = None,
                                 validate: Boolean = false): Map[String, ParameterLink] = {
    if (useManifests) {
      Map(
          Constants.OutputManifest -> ParameterLinkExec(execution,
                                                        Constants.OutputManifest,
                                                        Type.TFile)
      )
    } else {
      irOutputFields.map {
        case (fieldName, t) =>
          val fqn = prefix.map(p => s"${p}.${fieldName}").getOrElse(fieldName)
          if (validate) {
            validateOutput(Parameter.encodeName(fieldName), t)
          }
          fqn -> ParameterLinkExec(execution, fieldName, t)
      }
    }
  }

  def writeExecutionOutputLinks(execution: DxExecution, irOutputFields: Map[String, Type]): Unit = {
    val serializedLinks = createExecutionOutputLinks(execution, irOutputFields)
      .flatMap {
        case (name, link) => outputSerializer.createFieldsFromLink(link, name)
      }
      .filter {
        case (_, null | JsNull) => false
        case _                  => true
      }
    writeJsOutputs(serializedLinks)
  }

  def analysis: Option[DxAnalysis]

  def parentJob: Option[DxJob]

  def instanceType: Option[String]

  def getJobDetail(name: String): Option[JsValue]

  def getExecutableAttribute(name: String): Option[JsValue]

  def getExecutableDetail(name: String): Option[JsValue]

  lazy val language: Option[Language.Language] = getExecutableDetail(Constants.Language) match {
    case Some(JsString(lang)) => Some(Language.withName(lang))
    case None =>
      logger.warning("This applet ws built with an old version of dxCompiler - please rebuild")
      // we will attempt to detect the language/version later
      None
    case other =>
      throw new Exception(s"unexpected ${Constants.Language} value ${other}")
  }

  lazy val sourceCode: String = {
    val sourceCodeEncoded = getExecutableDetail(Constants.SourceCode) match {
      case Some(JsString(s)) => s
      case None =>
        logger.warning("This applet ws built with an old version of dxCompiler - please rebuild")
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

  lazy val targets: Option[Vector[String]] = {
    getExecutableDetail(Constants.Targets) match {
      case Some(JsString(t)) => Some(Vector(t))
      case Some(JsArray(a)) =>
        Some(a.map {
          case JsString(t) => t
          case other       => throw new Exception(s"invalid target value ${other}")
        })
      case None => None
      case other =>
        throw new Exception(s"invalid targets value ${other}")
    }
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
    getExecutableDetail(Constants.DelayWorkspaceDestruction) match {
      case Some(JsBoolean(flag)) => Some(flag)
      case None                  => None
    }

  lazy val blockPath: Vector[Int] = getExecutableDetail(Constants.BlockPath) match {
    case Some(JsArray(arr)) if arr.nonEmpty =>
      arr.map {
        case JsNumber(n) => n.toInt
        case _           => throw new Exception("Bad value ${arr}")
      }
    case Some(_: JsArray) | None => Vector.empty
    case other                   => throw new Exception(s"Bad value ${other}")
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

  lazy val useManifests: Boolean = getExecutableDetail(Constants.UseManifests) match {
    case Some(JsBoolean(b)) => b
    case None               => false
    case other =>
      throw new Exception(s"Invalid value ${other} for ${Constants.UseManifests}")
  }

  def error(e: Throwable): Unit
}

case class WorkerJobMeta(override val workerPaths: DxWorkerPaths = DxWorkerPaths.default,
                         override val dxApi: DxApi = DxApi.get,
                         override val logger: Logger = Logger.get)
    extends JobMeta(workerPaths, dxApi, logger) {
  lazy val project: DxProject = dxApi.currentProject

  private val rootDir = workerPaths.getRootDir()
  private val inputPath = rootDir.resolve(JobMeta.InputFile)
  private val outputPath = rootDir.resolve(JobMeta.OutputFile)
  private val jobInfoPath = rootDir.resolve(JobMeta.JobInfoFile)
  private val executableInfoPath = rootDir.resolve(JobMeta.ExecutableInfoFile)

  lazy override val rawJsInputs: Map[String, JsValue] = {
    if (Files.exists(inputPath)) {
      JsUtils.getFields(JsUtils.jsFromFile(inputPath))
    } else {
      logger.warning(s"input meta-file ${inputPath} does not exist")
      Map.empty
    }
  }

  def writeRawJsOutputs(outputJs: Map[String, JsValue]): Unit = {
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

  def jobId: String = dxApi.currentJob.id

  private lazy val jobDesc: DxJobDescribe = dxApi.currentJob.describe(
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
