package dx.executor

import java.nio.file.{Path, Paths}
import dx.api.{DxJob, InstanceTypeRequest}
import dx.core.getVersion
import dx.core.io.{
  DxdaManifest,
  DxdaManifestBuilder,
  DxfuseManifest,
  DxfuseManifestBuilder,
  StreamFiles
}
import dx.core.ir.{Type, TypeSerde, Value, ValueSerde}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.protocols.DxFileSource
import spray.json._
import dx.util.{
  AddressableFileNode,
  AddressableFileSource,
  Enum,
  FileUtils,
  LocalFileSource,
  SafeLocalizationDisambiguator,
  SysUtils,
  TraceLevel
}

object TaskAction extends Enum {
  type TaskAction = Value
  val CheckInstanceType, Prolog, InstantiateCommand, Epilog, Relaunch = Value
}

object TaskExecutor {
  val MaxDisambiguationDirs: Int = 5000
}

abstract class TaskExecutor(jobMeta: JobMeta,
                            streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
                            traceLengthLimit: Int = 10000) {

  private val fileResolver = jobMeta.fileResolver
  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger

  protected def trace(msg: String, minLevel: Int = TraceLevel.Verbose): Unit = {
    logger.traceLimited(msg, traceLengthLimit, minLevel)
  }

  private def printDirTree(): Unit = {
    if (logger.traceLevel >= TraceLevel.VVerbose) {
      trace("Directory structure:", TraceLevel.VVerbose)
      val (_, stdout, _) = SysUtils.execCommand("ls -lR", None)
      trace(stdout + "\n", TraceLevel.VVerbose)
    }
  }

  def executorName: String

  /**
    * Returns the minimal (i.e. cheapest) instance type that is
    * sufficient to the task's resource requirements.
    */
  protected def getInstanceTypeRequest: InstanceTypeRequest

  private def getRequestedInstanceType: String = {
    val instanceTypeRequest: InstanceTypeRequest = getInstanceTypeRequest
    logger.traceLimited(s"calcInstanceType $instanceTypeRequest")
    jobMeta.instanceTypeDb.apply(instanceTypeRequest).name
  }

  /**
    * Checks if we are already on the correct instance type. We only
    * relaunch if the current instance type is not sufficient for the
    * task requirements.
    */
  protected def checkInstanceType: Boolean = {
    // calculate the required instance type
    val requestedInstanceType = getRequestedInstanceType
    trace(s"requested instance type: ${requestedInstanceType}")
    val currentInstanceType = jobMeta.instanceType.getOrElse(
        throw new Exception(s"Cannot get instance type for job ${jobMeta.jobId}")
    )
    trace(s"current instance type: ${currentInstanceType}")
    val isSufficient =
      try {
        jobMeta.instanceTypeDb.matchesOrExceedes(currentInstanceType, requestedInstanceType)
      } catch {
        case ex: Throwable =>
          logger.warning("error comparing current and requested instance types",
                         exception = Some(ex))
          false
      }
    trace(s"isSufficient? ${isSufficient}")
    isSufficient
  }

  // These functions write/read the inputs and mappings between file
  // inputs and their local paths on disk. This is used to preserve
  // state between different phases of task execution, so we don't
  // have to re-evaluate expressions every time.

  private def writeEnv(schemas: Map[String, Type],
                       inputs: Map[String, (Type, Value)],
                       fileSourceToPath: Map[AddressableFileNode, Path]): Unit = {
    val schemasJs = schemas.values.foldLeft(Map.empty[String, JsValue]) {
      case (accu, schema) => TypeSerde.serialize(schema, accu)._2
    }
    val (inputsJs, newSchemasJs) =
      inputs
        .foldLeft((Map.empty[String, JsValue], schemasJs)) {
          case ((paramAccu, schemaAccu), (k, (t, v))) =>
            val (jsType, newSchemasJs) = TypeSerde.serialize(t, schemaAccu)
            val jsValue = ValueSerde.serialize(v)
            (paramAccu + (k -> JsObject("type" -> jsType, "value" -> jsValue)), newSchemasJs)
        }
    val uriToPath: Map[String, JsValue] = fileSourceToPath.map {
      case (fileSource: AddressableFileSource, path) =>
        fileSource.address -> JsString(path.toString)
      case (other, _) =>
        throw new RuntimeException(s"Can only serialize an AddressableFileSource, not ${other}")
    }
    val json = JsObject(
        "schemas" -> JsObject(newSchemasJs),
        "localizedInputs" -> JsObject(inputsJs),
        "dxUrlToPath" -> JsObject(uriToPath)
    )
    FileUtils.writeFileContent(jobMeta.workerPaths.getTaskEnvFile(), json.prettyPrint)
  }

  private def readEnv()
      : (Map[String, Type], Map[String, (Type, Value)], Map[AddressableFileNode, Path]) = {
    val (schemasJs, inputsJs, filesJs) =
      FileUtils.readFileContent(jobMeta.workerPaths.getTaskEnvFile()).parseJson match {
        case env: JsObject =>
          env.getFields("schemas", "localizedInputs", "dxUrlToPath") match {
            case Seq(JsObject(schemas), JsObject(inputs), JsObject(paths)) =>
              (schemas, inputs, paths)
            case _ =>
              throw new Exception("Malformed environment serialized to disk")
          }
        case _ => throw new Exception("Malformed environment serialized to disk")
      }
    val schemas = TypeSerde.deserializeSchemas(schemasJs)
    val (newSchemas, inputs) = inputsJs
      .foldLeft((schemas, Map.empty[String, (Type, Value)])) {
        case ((schemaAccu, paramAccu), (key, obj: JsObject)) =>
          obj.getFields("type", "value") match {
            case Seq(typeJs, valueJs) =>
              val (irType, newSchemas) = TypeSerde.deserialize(typeJs, schemaAccu)
              val irValue = ValueSerde.deserializeWithType(valueJs, irType, key)
              (newSchemas, paramAccu + (key -> (irType, irValue)))
            case _ =>
              throw new Exception("invalid env file")
          }
        case _ =>
          throw new Exception("invalid env file")
      }
    val fileSourceToPath = filesJs.map {
      case (uri, JsString(path)) => jobMeta.fileResolver.resolve(uri) -> Paths.get(path)
      case other                 => throw new Exception(s"unexpected path ${other}")
    }
    (newSchemas, inputs, fileSourceToPath)
  }

  /**
    * Returns the IR type and value for each task input, including default values
    * for any missing optional parameters.
    */
  protected def getInputsWithDefaults: Map[String, (Type, Value)]

  /**
    * Should we try to stream the file(s) associated with the given input parameter?
    */
  protected def streamFileForInput(parameterName: String): Boolean

  protected def getSchemas: Map[String, TSchema]

  def localizeInputFiles: (Map[String, (Type, Value)],
                           Map[AddressableFileNode, Path],
                           Option[DxdaManifest],
                           Option[DxfuseManifest]) = {

    assert(jobMeta.workerPaths.getInputFilesDir() != jobMeta.workerPaths.getDxfuseMountDir())

    val inputs = getInputsWithDefaults

    val (localFilesToPath, filesToStream, filesToDownload) =
      inputs.foldLeft(
          (Map.empty[AddressableFileNode, Path],
           Set.empty[AddressableFileNode],
           Set.empty[AddressableFileNode])
      ) {
        case ((localFilesToPath, filesToStream, filesToDownload), (name, (_, irValue))) =>
          val (local, remote) = extractFiles(irValue, fileResolver).foldLeft(
              (Map.empty[AddressableFileNode, Path], Set.empty[AddressableFileNode])
          ) {
            case ((local, remote), fs: LocalFileSource) =>
              // The file is already on the local disk, there is no need to download it.
              // TODO: make sure this file is NOT in the applet input/output directories.
              (local + (fs -> fs.canonicalPath), remote)
            case ((local, remote), other) =>
              (local, remote + other)
          }
          if (streamFiles == StreamFiles.All ||
              (streamFiles == StreamFiles.PerFile && streamFileForInput(name))) {
            (localFilesToPath ++ local, filesToStream ++ remote, filesToDownload)
          } else {
            (localFilesToPath ++ local, filesToStream, filesToDownload ++ remote)
          }
      }

    // build dxda and/or dxfuse manifests
    // We use a SafeLocalizationDisambiguator to determine the local path and deal
    // with file name collisions in the manner specified by the WDL spec. We set
    // separateDirsBySource = true because, when creating archives, we append one
    // directory at a time to the archive and then delete the source files, so the
    // smaller each append operation is, the less disk overhead is required.

    logger.traceLimited(s"downloading files = ${filesToDownload}")
    val downloadLocalizer =
      SafeLocalizationDisambiguator(
          jobMeta.workerPaths.getInputFilesDir(),
          existingPaths = localFilesToPath.values.toSet,
          separateDirsBySource = true,
          createDirs = true,
          disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
          logger = logger
      )
    val downloadFileSourceToPath: Map[AddressableFileNode, Path] =
      downloadLocalizer.getLocalPaths(filesToDownload)
    // write the manifest for dxda, if there are files to download
    val dxdaManifest = DxdaManifestBuilder(dxApi)
      .apply(downloadFileSourceToPath.collect {
        case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
      })

    logger.traceLimited(s"streaming files = ${filesToStream}")
    val streamingLocalizer =
      SafeLocalizationDisambiguator(
          jobMeta.workerPaths.getDxfuseMountDir(),
          existingPaths = localFilesToPath.values.toSet,
          separateDirsBySource = true,
          createDirs = false,
          disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
          logger = logger
      )
    val streamFileSourceToPath: Map[AddressableFileNode, Path] =
      streamingLocalizer.getLocalPaths(filesToStream)
    // write the manifest for dxfuse, if there are files to stream
    val dxfuseManifest = DxfuseManifestBuilder(dxApi)
      .apply(streamFileSourceToPath.collect {
        case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
      }, jobMeta.workerPaths)

    val fileSourceToPath = localFilesToPath ++ downloadFileSourceToPath ++ streamFileSourceToPath

    val uriToPath: Map[String, String] = fileSourceToPath.map {
      case (dxFs: DxFileSource, path)     => dxFs.address -> path.toString
      case (localFs: LocalFileSource, p2) => localFs.originalPath.toString -> p2.toString
      case other                          => throw new RuntimeException(s"unsupported file source ${other}")
    }

    // Replace the URIs with local file paths
    def pathTranslator(v: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      val uri = (t, v) match {
        case (_, VFile(uri))             => Some(uri)
        case (Some(TFile), VString(uri)) => Some(uri)
        case _                           => None
      }
      uri.map { u =>
        uriToPath.get(u) match {
          case Some(localPath)  => VFile(localPath)
          case None if optional => VNull
          case _ =>
            throw new Exception(s"Did not localize file ${u}")
        }
      }
    }

    val localizedInputs = inputs.view.mapValues {
      case (t, v) => (t, Value.transform(v, Some(t), pathTranslator))
    }.toMap

    (localizedInputs, fileSourceToPath, dxdaManifest, dxfuseManifest)
  }

  /**
    * For any File- and Directory-typed inputs for which a value is provided, materialize
    * those files on the local file system. This could be via direct download, or by producing
    * dxda and/or dxfuse manifests.
    *
    * Input files are represented as dx URLs (dx://proj-xxxx:file-yyyy::/A/B/C.txt) instead of
    * local files (/home/C.txt). Files may be referenced any number of times but are only
    * downloaded once.
    */
  def prolog(): Unit = {
    if (logger.isVerbose) {
      trace(s"Prolog debugLevel=${logger.traceLevel}")
      trace(s"dxCompiler version: ${getVersion}")
      printDirTree()
      trace(s"Task source code:\n${jobMeta.sourceCode}", traceLengthLimit)
    }

    val (localizedInputs, fileSourceToPath, dxdaManifest, dxfuseManifest) = localizeInputFiles

    dxdaManifest.foreach {
      case DxdaManifest(manifestJs) =>
        FileUtils.writeFileContent(jobMeta.workerPaths.getDxdaManifestFile(),
                                   manifestJs.prettyPrint)
    }

    dxfuseManifest.foreach {
      case DxfuseManifest(manifestJs) =>
        FileUtils.writeFileContent(jobMeta.workerPaths.getDxfuseManifestFile(),
                                   manifestJs.prettyPrint)
    }

    writeEnv(getSchemas, localizedInputs, fileSourceToPath)
  }

  /**
    * Generates and writes command script(s) to disk.
    * @param localizedInputs task inputs with localized files
    * @return localizedInputs updated with any additional (non-input) variables that
    *         may be required to evaluate the outputs.
    */
  protected def writeCommandScript(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)]

  def instantiateCommand(): Unit = {
    val (schemas, localizedInputs, fileSourceToPath) = readEnv()
    logger.traceLimited(s"InstantiateCommand, env = ${localizedInputs}")
    // evaluate the command block and write the command script
    val updatedInputs = writeCommandScript(localizedInputs)
    // write the updated env to disk
    writeEnv(schemas, updatedInputs, fileSourceToPath)
  }

  /**
    * Evaluates the outputs of the task. Returns mapping of output parameter
    * name to (type, value) and to (Set of tags, and Map of properties), where tags
    * and properites only apply to output files (or collections thereof).
    * @param localizedInputs the job inputs, with files localized to the worker
    */
  protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): (Map[String, (Type, Value)], Map[String, (Set[String], Map[String, String])])

  def epilog(): Unit = {
    if (logger.isVerbose) {
      trace(s"Epilog debugLevel=${logger.traceLevel}")
      printDirTree()
    }
    val (_, localizedInputs, fileSourceToPath) = readEnv()
    val (localizedOutputs, tagsAndProperties) = evaluateOutputs(localizedInputs)

    // Extract files from the outputs and build a map of all the string values
    // in the output values that might map to the same (absolute) local path.
    // Some of the outputs may be files that were inputs (in `fileSourceToPath`) -
    // these do not need to be re-uploaded. The `localPath`s will be the same but
    // the `originalPath`s may be different.
    val inputAddresses = fileSourceToPath.keySet.map(_.address)
    val inputPaths = fileSourceToPath.values.toSet
    val localOutputFileSources = localizedOutputs.map {
      case (name, (irType, irValue)) =>
        name -> extractOutputFiles(name, irValue, irType, fileResolver)
          .collect {
            case local: LocalFileSource
                if !(
                    inputAddresses.contains(local.address) ||
                      inputPaths.contains(local.canonicalPath)
                ) =>
              local
          }
    }

    val delocalizingValueToPath =
      localOutputFileSources.values.flatten.map(local => local.address -> local.canonicalPath).toMap

    // upload the files, and map their local paths to their remote URIs
    val delocalizedPathToUri = jobMeta
      .uploadOutputFiles(localOutputFileSources.map {
        case (name, localFileSources) => name -> localFileSources.map(_.canonicalPath)
      }, tagsAndProperties)
      .map {
        case (path, dxFile) => path -> dxFile.asUri
      }

    // Replace the local paths in the output values with URIs. For files that
    // were inputs, we can resolve them using a mapping of input values to URIs;
    // for files that were generated on the worker, this requires two look-ups:
    // first to get the absoulte Path associated with the file value (which may
    // be relative or absolute), and second to get the URI associated with the
    // Path. Returns an Optional[String] because optional outputs may be null.
    val inputValueToUri = fileSourceToPath
      .collect {
        case (fs: AddressableFileNode, path) =>
          Map(fs.address -> fs.address, path.toString -> fs.address)
      }
      .flatten
      .toMap
    def resolveFileValue(value: String): Option[String] = {
      inputValueToUri
        .get(value)
        .orElse(delocalizingValueToPath.get(value) match {
          case Some(path) => delocalizedPathToUri.get(path)
          case _          => None
        })
    }

    def pathTranslator(v: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      val uri = (t, v) match {
        case (_, VFile(uri))             => Some(uri)
        case (Some(TFile), VString(uri)) => Some(uri)
        case _                           => None
      }
      uri.map { u =>
        resolveFileValue(u) match {
          case Some(uri)        => VFile(uri)
          case None if optional => VNull
          case None =>
            throw new Exception(s"Did not delocalize file ${u}")
        }
      }
    }

    val delocalizedOutputs = localizedOutputs.view.mapValues {
      case (t, v) => (t, Value.transform(v, Some(t), pathTranslator))
    }.toMap

    // serialize the outputs to the job output file
    jobMeta.writeOutputs(delocalizedOutputs)
  }

  /**
    * Returns a mapping of output field names IR types.
    */
  protected def outputTypes: Map[String, Type]

  /**
    * Launches a sub-job with the same inputs and the dynamically calculated
    * instance type.
    */
  def relaunch(): Unit = {
    // Run a sub-job with the "body" entry point, and the required instance type
    val dxSubJob: DxJob =
      jobMeta.dxApi.runSubJob("body",
                              Some(getRequestedInstanceType),
                              JsObject(jobMeta.rawJsInputs),
                              Vector.empty,
                              jobMeta.delayWorkspaceDestruction)
    jobMeta.writeExecutionOutputLinks(dxSubJob, outputTypes)
  }

  def apply(action: TaskAction.TaskAction): String = {
    try {
      // setup the utility directories that the task-runner employs
      jobMeta.workerPaths.createCleanDirs()

      if (action == TaskAction.CheckInstanceType) {
        // special operation to check if this task is on the right instance type
        checkInstanceType.toString
      } else {
        action match {
          case TaskAction.Prolog             => prolog()
          case TaskAction.InstantiateCommand => instantiateCommand()
          case TaskAction.Epilog             => epilog()
          case TaskAction.Relaunch           => relaunch()
          case _ =>
            throw new Exception(s"Invalid executor action ${action}")
        }
        s"success ${action}"
      }
    } catch {
      case e: Throwable =>
        jobMeta.error(e)
        throw e
    }
  }
}
