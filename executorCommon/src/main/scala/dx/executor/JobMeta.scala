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
  DxFindDataObjects,
  DxFindDataObjectsConstraints,
  DxJob,
  DxJobDescribe,
  DxPath,
  DxProject,
  DxProjectDescribe,
  DxState,
  Field,
  FileUpload,
  InstanceTypeDB
}
import dx.core.Constants
import dx.core.io.{DxWorkerPaths, DxdaManifestBuilder}
import dx.core.ir.RunSpec.InstanceType
import dx.core.ir.Value.VNull
import dx.core.ir.{
  DxName,
  DxNameFactory,
  ExecutableLink,
  ExecutableLinkDeserializer,
  Manifest,
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
import dx.executor.JobMeta.MaxManifestJsLength
import dx.util.{
  CodecUtils,
  FileSourceResolver,
  FileUtils,
  JsUtils,
  Logger,
  StdMode,
  SysUtils,
  TraceLevel
}
import dx.util.protocols.DxFileAccessProtocol
import spray.json._

import scala.annotation.tailrec

object JobMeta {
  val InputFile = "job_input.json"
  val OutputFile = "job_output.json"
  val ErrorFile = "job_error.json"
  // this file has all the information about the job that is available
  // at execution time
  val JobInfoFile = "dnanexus-job.json"
  // this file has all the information about the executable
  val ExecutableInfoFile = "dnanexus-executable.json"
  val MaxConcurrentUploads = 8
  // max subjob input manifest document size for which a hash is used
  val MaxManifestJsLength = 100_000
  // functions used to download manifest files
  val ManifestDownloadDxda = "manifest_download_dxda"
  val WorkflowManifestDownloadDxda = "workflow_manifest_download_dxda"

  /**
    * Report an error, since this is called from a bash script, we can't simply
    * raise an exception. Instead, we write the error to a standard JSON file (job_error.json).
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
abstract class JobMeta(val workerPaths: DxWorkerPaths,
                       val dxNameFactory: DxNameFactory,
                       val dxApi: DxApi,
                       val logger: Logger) {
  def project: DxProject

  lazy val projectDesc: DxProjectDescribe = project.describe()

  def jobId: String

  def runJobScriptFunction(name: String,
                           successCodes: Option[Set[Int]] = Some(Set(0)),
                           truncateLogs: Boolean = true,
                           forwardStd: Boolean = false): Unit

  def rawJsInputs: Map[DxName, JsValue]

  val manifestRootFolder = "/.d"
  lazy val manifestFolder = s"${manifestRootFolder}/${jobId}"
  lazy val manifestProjectAndFolder: String = s"${dxApi.currentProjectId.get}:${manifestFolder}"
  lazy val projectOutputFolder: String = {
    if (useManifests) {
      manifestFolder
    } else {
      folder
    }
  }

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
  private def unpackManifests(rawInputs: Map[DxName, JsValue]): Map[DxName, JsValue] = {
    logger.traceLimited(s"""Unpacking manifests from raw inputs:
                           |${JsObject(rawInputs.map {
                             case (dxName, jsv) => dxName.decoded -> jsv
                           }).prettyPrint}""".stripMargin)

    // resolve a file input to a DxFile
    def resolveManifestFiles(key: DxName): Vector[DxFile] = {
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
    val manifestLinks = (rawInputs.get(Constants.InputLinks) match {
      case Some(JsObject(fields)) if fields.contains(Constants.ComplexValueKey) =>
        fields(Constants.ComplexValueKey).asJsObject.fields
      case Some(JsObject(fields)) => fields
      case None                   => Map.empty[String, JsValue]
      case other =>
        throw new Exception(s"invalid manifest links ${other}")
    }).map {
      case (name, jsv) => dxNameFactory.fromDecodedName(name) -> jsv
    }

    if (manifestLinks.isEmpty) {
      // if there are no manifest links, then there must be at most one manifest
      // (either as a hash or a file). We can just parse it and return the values.
      val manifestJson = (manifestHash, manifestFiles) match {
        case (None, Vector())         => None
        case (Some(inputs), Vector()) => Some(inputs)
        case (None, Vector(manifestFile)) =>
          Some(new String(dxApi.downloadBytes(manifestFile, retryLimit = 10)).parseJson)
        case _ =>
          throw new Exception(
              "manifest links are required when there is more than one manifest file"
          )
      }
      manifestJson
        .map(Manifest.parse(_, dxNameFactory).jsValues)
        .getOrElse(Map.empty[DxName, JsValue])
    } else {
      def downloadAndParseManifests(files: Vector[DxFile],
                                    dxdaManifestFile: Path,
                                    downloadFunction: String): Vector[Manifest] = {
        if (files.isEmpty) {
          Vector.empty
        } else if (files.size == 1) {
          // download manifest bytes, convert to JSON, and parse into Manifest object
          Vector(
              Manifest.parse(new String(dxApi.downloadBytes(files.head, retryLimit = 10)).parseJson,
                             dxNameFactory)
          )
        } else {
          // bulk describe manifest files
          val manifestFileToPath = DxFindDataObjects(dxApi, logger = logger)
            .query(
                DxFindDataObjectsConstraints(
                    project = dxApi.currentProject,
                    folder = Some(manifestRootFolder),
                    recurse = true,
                    objectClass = Some("file"),
                    ids = files.map(_.id).toSet
                ),
                describe = true,
                defaultFields = true,
                extraFields = Set(Field.Parts)
            )
            .keys
            .map {
              case dxFile: DxFile =>
                dxFile -> workerPaths.getManifestFilesDir().resolve(dxFile.getName).asJavaPath
              case other => throw new Exception(s"unexpected result ${other}")
            }
            .toMap
          if (manifestFileToPath.size != files.size) {
            throw new Exception(
                s"""failed to describe one or more manifest file(s)
                   |${files}
                   |vs
                   |${manifestFileToPath}""".stripMargin
            )
          }
          // write the dxda manifest to a file
          FileUtils.writeFileContent(
              dxdaManifestFile,
              DxdaManifestBuilder(dxApi, logger).apply(manifestFileToPath).get.value.prettyPrint
          )
          // run dxda via a subprocess
          runJobScriptFunction(downloadFunction)
          // read the manifest files
          manifestFileToPath.values.map { localManifestFile =>
            val manifest = Manifest.parseFile(localManifestFile, dxNameFactory)
            // there is more than one manifest so check that they all have IDs
            if (manifest.id.isEmpty) {
              throw new Exception("when there are multiple manifests, all manifests must have ID")
            }
            manifest
          }.toVector
        }
      }

      // parse all the manifests and combine into a single Vector
      val manifests = manifestHash
        .map(h => Vector(Manifest.parse(h, dxNameFactory)))
        .getOrElse(Vector.empty) ++
        downloadAndParseManifests(manifestFiles,
                                  workerPaths.getDxdaManifestDownloadManifestFile().asJavaPath,
                                  JobMeta.ManifestDownloadDxda)
      // create lookup function for manifests
      val manifestLookup: (String, DxName) => Option[JsValue] = manifests.size match {
        case 0 => (_: String, _: DxName) => None
        case 1 =>
          val values = manifests.head.jsValues
          (_: String, manifestParamName: DxName) => values.get(manifestParamName)
        case _ =>
          val manifestMap = manifests.map {
            case m: Manifest if m.id.isDefined => m.id.get -> m.jsValues
            case _ =>
              throw new Exception(
                  "when there are multiple manifests, all manifests must have id defined"
              )
          }.toMap
          (manifestId: String, manifestParamName: DxName) => {
            manifestMap.get(manifestId).flatMap(_.get(manifestParamName))
          }
      }

      // create lookup function for workflow manifests - this is lazy because the
      // manifest links may not contain any workflow references
      lazy val workflowManifestLookup: DxName => Option[JsValue] = {
        val workflowManifestHash = rawInputs.get(Constants.WorkflowInputManifest) match {
          case Some(hash: JsObject) => Some(hash)
          case None                 => None
          case other =>
            throw new Exception(s"invalid manifest ${other}")
        }
        val workflowManifestFiles = resolveManifestFiles(Constants.WorkflowInputManifestFiles)
        val workflowManifests = workflowManifestHash
          .map(h => Vector(Manifest.parse(h, dxNameFactory)))
          .getOrElse(Vector.empty) ++
          downloadAndParseManifests(
              workflowManifestFiles,
              workerPaths.getDxdaWorkflowManifestDownloadManifestFile().asJavaPath,
              JobMeta.WorkflowManifestDownloadDxda
          )
        val workflowManifestLinks: Map[DxName, JsValue] =
          (rawInputs.get(Constants.WorkflowInputLinks) match {
            case Some(JsObject(fields)) if fields.contains(Constants.ComplexValueKey) =>
              fields(Constants.ComplexValueKey).asJsObject.fields
            case Some(JsObject(fields)) => fields
            case None                   => Map.empty[String, JsValue]
            case other =>
              throw new Exception(s"invalid workflow manifest links ${other}")
          }).map {
            case (name, jsv) => dxNameFactory.fromDecodedName(name) -> jsv
          }
        if (workflowManifests.isEmpty) {
          throw new Exception("there are no workflow manifest files")
        } else if (workflowManifestLinks.isEmpty) {
          if (workflowManifests.size == 1) {
            val values = workflowManifests.head.jsValues
            (paramName: DxName) => values.get(paramName)
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
          def workflowManifestLookup(paramName: DxName): Option[JsValue] = {
            workflowManifestLinks.get(paramName) match {
              case Some(JsObject(fields)) if fields.keySet == Set(ValueSerde.WrappedValueKey) =>
                fields(ValueSerde.WrappedValueKey).asJsObject.fields.head match {
                  case (key, value) if key == Constants.ValueKey.decoded => Some(value)
                  case (key, JsString(workflowParamName)) if key == Constants.WorkflowKey.decoded =>
                    workflowManifestLookup(dxNameFactory.fromDecodedName(workflowParamName))
                  case (manifestId, JsString(manifestParamName)) =>
                    workflowManifestMap
                      .get(manifestId)
                      .flatMap(_.get(dxNameFactory.fromDecodedName(manifestParamName)))
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
            case (key, value) if key == Constants.ValueKey.decoded => Some(value)
            case (key, JsString(workflowParamName)) if key == Constants.WorkflowKey.decoded =>
              workflowManifestLookup(dxNameFactory.fromDecodedName(workflowParamName))
            case (manifestId, JsString(manifestParamName)) =>
              manifestLookup(manifestId, dxNameFactory.fromDecodedName(manifestParamName))
            case other =>
              throw new Exception(s"invalid manifest link entry ${other}")
          }
          value.map(paramName -> _)
        case other =>
          throw new Exception(s"invalid manifest link ${other}")
      }
    }
  }

  lazy val (jsInputs: Map[DxName, JsValue], jsOverrides: Option[JsValue]) = {
    // pop the overrides off the rest of the inputs
    val (rawOverrides, rawInputs) = rawJsInputs.partition {
      case (dxName, _) => dxName == Constants.Overrides
    }
    val jsInputs = if (useManifests) {
      unpackManifests(rawInputs)
    } else {
      rawInputs
    }
    val jsOverrides = rawOverrides.values.headOption.map {
      case JsObject(fields) if fields.contains(Constants.ComplexValueKey) =>
        fields(Constants.ComplexValueKey)
      case jsv => jsv
    }
    if (logger.isVerbose) {
      if (jsInputs.isEmpty) {
        logger.trace("No inputs")
      } else {
        logger.traceLimited(s"Raw inputs:\n${JsObject(jsInputs.map {
          case (dxName, jsv) => dxName.decoded -> jsv
        }).prettyPrint}")
      }
      if (jsOverrides.isEmpty) {
        logger.trace("No overrides")
      } else {
        logger.traceLimited(s"Raw overrides:\n${jsOverrides.get.prettyPrint}")
      }
    }
    (jsInputs, jsOverrides)
  }

  private lazy val allFilesReferenced: Vector[DxFile] = {
    // bulk describe all files referenced in the inputs
    logger.trace("Discovering all files in the input values")
    val queryFiles = jsInputs.values.flatMap(DxFile.findFiles(dxApi, _)).toVector.distinct
    logger.trace(s"Bulk describing ${queryFiles.size} files")
    val dxFiles = dxApi.describeFilesBulk(queryFiles, searchWorkspaceFirst = true, validate = true)
    // check that all files are in the closed state
    logger.trace(s"Checking that all files are closed")
    val notClosed = dxFiles.filterNot(_.describe().state == DxState.Closed)
    if (notClosed.nonEmpty) {
      throw new Exception(
          s"input file(s) not in the 'closed' state: ${notClosed.map(_.id).mkString(",")}"
      )
    }
    logger.trace(s"Successfully described ${dxFiles.size} files")
    dxFiles
  }

  protected lazy val dxFileDescCache: DxFileDescCache = DxFileDescCache(allFilesReferenced)

  lazy val inputDeserializer: ParameterLinkDeserializer =
    ParameterLinkDeserializer(dxFileDescCache, dxApi)

  lazy val inputSpec: Map[String, Type] = {
    val inputSpec: Map[String, Type] = getExecutableAttribute("inputSpec")
      .map {
        case JsArray(spec) => TypeSerde.fromNativeSpec(spec)
        case other         => throw new Exception(s"invalid inputSpec ${other}")
      }
      .getOrElse(Map.empty)
    if (logger.isVerbose) {
      logger.traceLimited(s"inputSpec:\n  ${inputSpec.mkString("\n  ")}")
    }
    inputSpec
  }

  lazy val inputs: Map[DxName, Value] = {
    // If we are using manifests, then the inputSpec won't match the task inputs.
    // Otherwise, if we have access to the inputSpec, use it to guide deserialization.
    val inputs = jsInputs.map {
      case (dxName, value) if useManifests =>
        dxName -> inputDeserializer.deserializeInput(value)
      case (dxName, value) =>
        val irValue = inputSpec.get(dxName.encoded) match {
          case None =>
            logger.warning(s"inputSpec is missing field ${dxName}")
            inputDeserializer.deserializeInput(value)
          case Some(t) if Type.unwrapOptional(t) == Type.THash =>
            logger.trace(
                s"""expected type of input field '${dxName}' is THash, which may represent an
                   |unknown schema type, so deserializing without type""".stripMargin
                  .replaceAll("\n", " ")
            )
            inputDeserializer.deserializeInput(value)
          case Some(t) =>
            inputDeserializer.deserializeInputWithType(value, t, dxName.decoded)
        }
        dxName -> irValue
    }
    if (logger.isVerbose && inputs.nonEmpty) {
      val inputStr = inputs
        .map {
          case (name, value) => s"${name}: ${ValueSerde.toString(value)}"
        }
        .mkString("\n  ")
      logger.traceLimited(s"Deserialized inputs:\n  ${inputStr}")
    }
    inputs
  }

  lazy val primaryInputs: Map[DxName, Value] = {
    // discard auxiliary fields
    inputs.filter {
      case (dxName, _) => !dxName.suffix.exists(_.endsWith(Constants.FlatFilesSuffix))
    }
  }

  lazy val fileResolver: FileSourceResolver = {
    logger.trace(s"Creating FileSourceResolver localDirectories = ${workerPaths.getWorkDir()}")
    val dxProtocol = DxFileAccessProtocol(dxApi, dxFileDescCache)
    val fileResolver = FileSourceResolver.create(
        localDirectories = Vector(workerPaths.getWorkDir().asJavaPath),
        userProtocols = Vector(dxProtocol),
        logger = logger
    )
    FileSourceResolver.set(fileResolver)
    fileResolver
  }

  def uploadFiles(filesToUpload: Iterable[FileUpload]): Map[Path, DxFile]

  protected def writeRawJsOutputs(outputJs: Map[DxName, JsValue]): Unit

  def writeJsOutputs(outputJs: Map[DxName, JsValue],
                     skipExtraManifestOutputs: Boolean = false): Unit = {
    val rawOutputJs = if (useManifests && !outputJs.contains(Constants.OutputManifest)) {
      val manifestId = rawJsInputs.get(Constants.OutputId) match {
        case Some(JsString(id)) => id
        case other =>
          throw new Exception(s"missing or invalid outputId ${other}")
      }
      val manifestValues = rawJsInputs.get(Constants.CallName) match {
        case Some(JsString(callName)) =>
          outputJs.map {
            case (dxName, value) => dxName.pushDecodedNamespace(callName) -> value
          }
        case None => outputJs
        case other =>
          throw new Exception(s"invalid callName ${other}")
      }
      val extraManifestOutputJs = Option
        .when(!skipExtraManifestOutputs) {
          rawJsInputs.get(Constants.ExtraOutputs).map {
            case JsObject(fields) =>
              fields.map {
                case (key, value) => dxNameFactory.fromDecodedName(key) -> value
              }
            case other => throw new Exception(s"invalid extra_outputs___ value ${other}")
          }
        }
        .flatten
        .getOrElse(Map.empty)
      val manifest = Manifest(manifestValues ++ extraManifestOutputJs, id = Some(manifestId))
      val destination = s"${manifestProjectAndFolder}/${jobId}_output.manifest.json"
      val manifestDxFile = dxApi.uploadString(manifest.toJson().prettyPrint, destination)
      outputSerializer
        .createFields(Constants.OutputManifest, Type.TFile, Value.VFile(manifestDxFile.asUri))
        .toMap
    } else {
      outputJs
    }
    writeRawJsOutputs(rawOutputJs)
  }

  lazy val outputSerializer: ParameterLinkSerializer =
    ParameterLinkSerializer(fileResolver, dxApi, pathsAsObjects = pathsAsObjects)

  lazy val executableLinkDeserializer: ExecutableLinkDeserializer =
    ExecutableLinkDeserializer(dxNameFactory, dxApi)

  private val detailToFilenameRegex = "[^\\w_]".r

  /**
    * Prepares the inputs for a subjob. If using manifests, creates a manifest
    * with the same output ID as the current job.
    * @param inputs raw subjob inputs
    * @param extraManifestOutputs extra outputs that need to be added to the output manifest of the
    *                             called job; ignored if useManifests=false
    * @param executableLink the subjob executable
    * @param callName the (unencoded) call name to use to prefix the subjob outputs
    * @return serialized subjob inputs
    */
  def prepareSubjobInputs(inputs: Map[DxName, (Type, Value)],
                          extraManifestOutputs: Option[Map[DxName, (Type, Value)]] = None,
                          executableLink: ExecutableLink,
                          callName: Option[String] = None,
                          nameDetail: Option[String] = None): Map[DxName, JsValue] = {
    val inputsJs = outputSerializer.createFieldsFromMap(inputs)
    // Check that we have all the compulsory arguments.
    // Note that we don't have the information here to tell difference between optional and non-
    // optional parameters, so we emit warnings for missing arguments.
    executableLink.inputs.keys.filterNot(inputsJs.contains).foreach { argName =>
      logger.warning(s"Missing argument ${argName} to call ${executableLink.name}", force = true)
    }
    if (useManifests) {
      val manifestValuesJs = JsObject(inputsJs.map {
        case (dxName, jsv) => dxName.decoded -> jsv
      })
      val manifestJsStr = manifestValuesJs.compactPrint
      val manifestInputJs = if (manifestJsStr.length <= MaxManifestJsLength) {
        Map(Constants.InputManifest -> manifestValuesJs)
      } else {
        // For large subjob inputs, put them in a compact manifest file and upload it.
        // If there is a name detail, there are probably going to be a large number of manifest
        // files, so put them in a subfolder.
        // TODO: it is unclear if we'll ever need a manifest id here, and if so, what the id should be
        val manifestDetail =
          nameDetail.map(d => s"/${detailToFilenameRegex.replaceAllIn(d, "_")}").getOrElse("")
        val manifestFilename =
          s"${manifestProjectAndFolder}/${jobId}_subjob${manifestDetail}.manifest.json"
        val manifestDxFile =
          dxApi.uploadString(manifestJsStr, manifestFilename)
        Map(Constants.InputManifestFiles -> JsArray(manifestDxFile.asJson))
      }
      val commonInputsJs = Vector(
          Some(Constants.OutputId -> rawJsInputs(Constants.OutputId)),
          extraManifestOutputs.map { extraOutputs =>
            Constants.ExtraOutputs -> JsObject(
                outputSerializer.createFieldsFromMap(extraOutputs).map {
                  case (dxName, jsv) => dxName.decoded -> jsv
                }
            )
          },
          callName.map(name => Constants.CallName -> JsString(name))
      ).flatten.toMap
      manifestInputJs ++ commonInputsJs
    } else {
      inputsJs
    }
  }

  @tailrec
  private def outputTypesEqual(expectedType: Type, actualType: Type): Boolean = {
    (expectedType, actualType) match {
      case (a, b) if a == b             => true
      case (a, b) if Type.isOptional(a) =>
        // non-optional actual type is compatible with optional expected type
        outputTypesEqual(Type.unwrapOptional(a), Type.unwrapOptional(b))
      case (Type.TArray(a, false), Type.TArray(b, true)) =>
        // non-empty actual array is compatible with maybe-empty expected array
        outputTypesEqual(a, b)
      case (a, b) if !Type.isNative(b) =>
        // non-native types are represented as hashes
        a == Type.THash
      case _ => false
    }
  }

  private lazy val outputSpec: Map[String, Type] = {
    val outputSpec: Map[String, Type] = getExecutableAttribute("outputSpec")
      .map {
        case JsArray(spec) => TypeSerde.fromNativeSpec(spec)
        case other         => throw new Exception(s"invalid outputSpec ${other}")
      }
      .getOrElse(Map.empty)
    logger.traceLimited(s"outputSpec:\n  ${outputSpec.mkString("\n  ")}")
    outputSpec
  }

  private def validateOutput(dxName: DxName,
                             actualType: Type,
                             value: Option[Value] = None): Type = {
    outputSpec.get(dxName.encoded) match {
      case Some(expectedType) if outputTypesEqual(expectedType, actualType) =>
        actualType
      case Some(expectedType) if Type.unwrapOptional(expectedType) == Type.THash =>
        logger.trace(
            s"""Expected type of output field ${dxName.encoded} is THash which may represent an
               |unknown schema type, so deserializing without type""".stripMargin
              .replaceAll("\n", " ")
        )
        actualType
      case Some(expectedType) if value.isDefined =>
        logger.trace(s"Coercing value ${value.get} to ${expectedType}")
        try {
          // TODO: add a Type.isCoercibleTo function
          logger.ignore(Value.coerceTo(value.get, expectedType))
          expectedType
        } catch {
          case _: Throwable =>
            throw new Exception(
                s"""value ${value} is not coercible to expected type ${expectedType}""".stripMargin
                  .replaceAll("\n", " ")
            )
        }
      case Some(expectedType) =>
        throw new Exception(
            s"""output field ${dxName.encoded} has mismatch between actual type ${actualType}
               |and expected type ${expectedType}""".stripMargin.replaceAll("\n", " ")
        )
      case None =>
        logger.warning(s"outputSpec is missing field ${dxName.encoded}")
        actualType
    }
  }

  def writeOutputs(outputs: Map[DxName, (Type, Value)]): Unit = {
    // If we are using manifests, then the outputSpec won't match the task outputs.
    // Otherwise, if we have access to the outputSpec, use it to validate the outputs.
    if (!useManifests) {
      outputs.foreach {
        case (name, (actualType, value)) => validateOutput(name, actualType, Some(value))
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

  def createOutputLinks(outputs: Map[DxName, (Type, Value)],
                        validate: Boolean = true): Map[DxName, ParameterLink] = {
    outputs.collect {
      case (name, (actualType, value)) if value != VNull =>
        val validatedType = if (validate) {
          validateOutput(name, actualType, Some(value))
        } else {
          actualType
        }
        name -> outputSerializer.createLink(validatedType, value)
    }
  }

  def writeOutputLinks(outputs: Map[DxName, ParameterLink]): Unit = {
    val serializedLinks = outputs
      .map {
        case (name, link) => outputSerializer.createFieldsFromLink(link, name)
      }
      .flatten
      .toMap
    writeJsOutputs(serializedLinks)
  }

  def createExecutionOutputLinks(execution: DxExecution,
                                 irOutputFields: Map[DxName, Type],
                                 prefix: Option[String] = None,
                                 validate: Boolean = false): Map[DxName, ParameterLink] = {
    if (useManifests) {
      Map(
          Constants.OutputManifest -> ParameterLinkExec.create(execution,
                                                               Constants.OutputManifest,
                                                               Type.TFile)
      )
    } else {
      irOutputFields.map {
        case (fieldDxName, actualType) =>
          val fqn = prefix
            .map(fieldDxName.pushDecodedNamespace)
            .getOrElse(fieldDxName)
          val validatedType = if (validate) {
            validateOutput(fieldDxName, actualType)
          } else {
            actualType
          }
          fqn -> ParameterLinkExec.create(execution, fieldDxName, validatedType)
      }
    }
  }

  def writeExecutionOutputLinks(execution: DxExecution, irOutputFields: Map[DxName, Type]): Unit = {
    val serializedLinks = createExecutionOutputLinks(execution, irOutputFields)
      .flatMap {
        case (name, link) => outputSerializer.createFieldsFromLink(link, name)
      }
      .filter {
        case (_, null | JsNull) => false
        case _                  => true
      }
    writeJsOutputs(serializedLinks, skipExtraManifestOutputs = true)
  }

  /**
    * The parent analysis of the current job.
    */
  def analysis: Option[DxAnalysis]

  /**
    * The parent job of the current job.
    */
  def parentJob: Option[DxJob]

  /**
    * The current instance type.
    */
  def instanceType: Option[String]

  /**
    * The job output folder.
    */
  def folder: String

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

  lazy val parserOptions: Option[JsValue] = getExecutableDetail(Constants.ParseOptions)

  lazy val instanceTypeDb: InstanceTypeDB = {
    try {
      // first try to create the db from the current project
      InstanceType.createDb(dxApi = dxApi)
    } catch {
      case _: Throwable =>
        // fall back to the cached database if available
        getExecutableDetail(Constants.InstanceTypeDb) match {
          case Some(JsString(s)) =>
            val js = CodecUtils.base64DecodeAndGunzip(s)
            js.parseJson.convertTo[InstanceTypeDB]
          case other =>
            throw new Exception(s"unexpected ${Constants.InstanceTypeDb} value ${other}")
        }
    }
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
      case other =>
        throw new Exception(s"invalid ${Constants.DelayWorkspaceDestruction} value ${other}")
    }

  lazy val blockPath: Vector[Int] = getExecutableDetail(Constants.BlockPath) match {
    case Some(JsArray(arr)) if arr.nonEmpty =>
      arr.map {
        case JsNumber(n) => n.toInt
        case other =>
          throw new Exception(s"Invalid array item ${other} for ${Constants.BlockPath}")
      }
    case Some(_: JsArray) | None => Vector.empty
    case other =>
      throw new Exception(s"Invalid value ${other} for ${Constants.BlockPath}")
  }

  lazy val scatterStart: Int = getJobDetail(Constants.ContinueStart) match {
    case Some(JsNumber(s)) => s.toIntExact
    case Some(other) =>
      throw new Exception(s"Invalid value ${other} for ${Constants.ContinueStart}")
    case _ => 0
  }

  lazy val scatterSize: Int = getExecutableDetail(Constants.ScatterChunkSize) match {
    case Some(JsNumber(n)) => n.toIntExact
    case None              => Constants.JobPerScatterDefault
    case other =>
      throw new Exception(s"Invalid value ${other} for ${Constants.ScatterChunkSize}")
  }

  lazy val scatterOutputShape: Option[Vector[Int]] = getJobDetail(Constants.OutputShape).map {
    case JsArray(arr) =>
      arr.map {
        case JsNumber(i) => i.toIntExact
        case other       => throw new Exception(s"invalid output shape value ${other}")
      }
    case other =>
      throw new Exception(s"invalid ${Constants.OutputShape} value ${other}")
  }

  lazy val scatterSkippedIndices: Option[Set[Int]] = getJobDetail(Constants.SkippedIndices).map {
    case JsArray(indices) =>
      indices.map {
        case JsNumber(n) if n.isValidInt => n.toIntExact
        case other =>
          throw new Exception(s"invalid skipped index value ${other}")
      }.toSet
    case other =>
      throw new Exception(s"invalid ${Constants.SkippedIndices} value ${other}")
  }

  lazy val useManifests: Boolean = getExecutableDetail(Constants.UseManifests) match {
    case Some(JsBoolean(b)) => b
    case None               => false
    case other =>
      throw new Exception(s"Invalid value ${other} for ${Constants.UseManifests}")
  }

  lazy val pathsAsObjects: Boolean = getExecutableDetail(Constants.PathsAsObjects) match {
    case Some(JsBoolean(b)) => b
    case None               => false
    case other =>
      throw new Exception(s"Invalid value ${other} for ${Constants.PathsAsObjects}")
  }

  lazy val fileDependencies: Set[String] = {
    getExecutableDetail(Constants.FileDependencies) match {
      case Some(JsArray(deps)) =>
        deps.map {
          case JsString(uri) => uri
          case other =>
            throw new Exception(s"invalid ${Constants.FileDependencies} value: ${other}")
        }.toSet
      case None => Set.empty
      case other =>
        throw new Exception(s"invalid value for ${Constants.FileDependencies}: ${other}")
    }
  }

  def error(e: Throwable): Unit
}

case class WorkerJobMeta(override val workerPaths: DxWorkerPaths,
                         waitOnUpload: Boolean,
                         override val dxNameFactory: DxNameFactory,
                         override val dxApi: DxApi = DxApi.get,
                         override val logger: Logger = Logger.get)
    extends JobMeta(workerPaths, dxNameFactory, dxApi, logger) {
  lazy val project: DxProject =
    dxApi.currentProject.getOrElse(throw new Exception("no current project"))

  private val rootDir = workerPaths.getRootDir().asJavaPath
  private val inputPath = rootDir.resolve(JobMeta.InputFile)
  private val outputPath = rootDir.resolve(JobMeta.OutputFile)
  private val jobInfoPath = rootDir.resolve(JobMeta.JobInfoFile)
  private val executableInfoPath = rootDir.resolve(JobMeta.ExecutableInfoFile)
  private val LogLimit = 1_000_000
  private val stdLimit = Map("lines" -> 10, "chars" -> 500)

  def codeFile: Path = {
    workerPaths.getRootDir().resolve(s"${jobId}.code.sh").asJavaPath
  }

  override def runJobScriptFunction(name: String,
                                    successCodes: Option[Set[Int]] = Some(Set(0)),
                                    truncateLogs: Boolean = true,
                                    forwardStd: Boolean = false): Unit = {
    val command = s"bash -c 'source ${codeFile} && ${name}'"
    logger.trace(s"Running job script function ${name}")
    if (forwardStd) {
      val (rc, _, _) = SysUtils.runCommand(command,
                                           exceptionOnFailure = false,
                                           stdoutMode = StdMode.Forward,
                                           stderrMode = StdMode.Forward)
      if (!successCodes.forall(_.contains(rc))) {
        val subprocessStderr = FileUtils.readFileContent(workerPaths.getStderrFile().asJavaPath)
        throw new Exception(s"""job script function ${name} exited with permanent fail code ${rc}
                               |${truncateStd(subprocessStderr)}
                               |""".stripMargin)
      }
    } else {
      val (rc, stdout, stderr) = SysUtils.execCommand(command, exceptionOnFailure = false)
      if (successCodes.forall(_.contains(rc))) {
        val limit = if (truncateLogs) Some(LogLimit) else None
        logger.trace(
            s"""Job script function ${name} exited with success code ${rc}
               |----- stdout -----
               |${stdout}
               |------------------""".stripMargin,
            maxLength = limit,
            showBeginning = true,
            showEnd = true
        )
      } else {
        logger.error(s"""Job script function ${name} exited with permanent fail code ${rc}
                        |----- stdout -----:
                        |${stdout}
                        |----- stderr-----:
                        |${stderr}
                        |-----------------""".stripMargin)
        throw new Exception(s"job script function ${name} exited with permanent fail code ${rc}")
      }
    }
  }

  private def truncateStd(stdString: String): String = {
    val lastLines = stdString
      .split("\n")
      .takeRight(stdLimit("lines"))
      .mkString("\n")
    lastLines.takeRight(Math.min(lastLines.length(), stdLimit("chars")))

  }

  lazy override val rawJsInputs: Map[DxName, JsValue] = {
    if (Files.exists(inputPath)) {
      logger.trace(s"Loading raw JSON input from ${inputPath}")
      JsUtils.getFields(JsUtils.jsFromFile(inputPath)).map {
        case (name, jsv) => dxNameFactory.fromEncodedName(name) -> jsv
      }
    } else {
      logger.warning(s"input meta-file ${inputPath} does not exist")
      Map.empty
    }
  }

  override def uploadFiles(filesToUpload: Iterable[FileUpload]): Map[Path, DxFile] = {
    dxApi.uploadFiles(filesToUpload, waitOnUpload, maxConcurrent = JobMeta.MaxConcurrentUploads)
  }

  def writeRawJsOutputs(outputJs: Map[DxName, JsValue]): Unit = {
    JsUtils.jsToFile(JsObject(outputJs.map {
      case (dxName, jsv) => dxName.encoded -> jsv
    }), outputPath)
  }

  private lazy val jobInfo: Map[String, JsValue] = {
    if (Files.exists(jobInfoPath)) {
      FileUtils.readFileContent(jobInfoPath).parseJson.asJsObject.fields
    } else {
      logger.warning(s"info meta-file ${jobInfoPath} does not exist")
      Map.empty
    }
  }

  def jobId: String = dxApi.currentJobId.get

  private lazy val jobDesc: DxJobDescribe = dxApi.currentJob.get.describe(
      Set(Field.Details, Field.Executable, Field.Folder, Field.InstanceType)
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

  def folder: String = jobDesc.folder.get

  def getJobDetail(name: String): Option[JsValue] = jobDetails.get(name)

  def getExecutableDetail(name: String): Option[JsValue] = executableDetails.get(name)

  def error(e: Throwable): Unit = {
    JobMeta.writeError(rootDir, e)
  }
}
