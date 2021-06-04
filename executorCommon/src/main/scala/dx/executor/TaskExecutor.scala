package dx.executor

import java.nio.file.{Files, Path, Paths}
import dx.api.{DxJob, DxPath, InstanceTypeRequest}
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
import dx.util.protocols.{DxFileSource, DxFolderSource}
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
                            fileUploader: FileUploader = SerialFileUploader(),
                            streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
                            waitOnUpload: Boolean = false,
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
                       fileSourceToPath: Map[AddressableFileNode, Path],
                       folderSourceToPath: Map[AddressableFileSource, Path]): Unit = {
    val schemasJs = schemas.values.foldLeft(Map.empty[String, JsValue]) {
      case (accu, schema) => TypeSerde.serialize(schema, accu)._2
    }
    val (inputsJs, newSchemasJs) =
      inputs
        .foldLeft((Map.empty[String, JsValue], schemasJs)) {
          case ((paramAccu, schemaAccu), (k, (t, v))) =>
            val (jsType, newSchemasJs) = TypeSerde.serialize(t, schemaAccu)
            val jsValue = ValueSerde.serialize(v, pathsAsObjects = jobMeta.pathsAsObjects)
            (paramAccu + (k -> JsObject("type" -> jsType, "value" -> jsValue)), newSchemasJs)
        }
    val fileUriToPath: Map[String, JsValue] = fileSourceToPath.map {
      case (fileSource: AddressableFileSource, path) =>
        fileSource.address -> JsString(path.toString)
      case (other, _) =>
        throw new RuntimeException(s"Can only serialize an AddressableFileSource, not ${other}")
    }
    val folderUriToPath: Map[String, JsValue] = folderSourceToPath.map {
      case (folderSource: AddressableFileSource, path) =>
        folderSource.address -> JsString(path.toString)
    }
    val json = JsObject(
        "schemas" -> JsObject(newSchemasJs),
        "localizedInputs" -> JsObject(inputsJs),
        "fileDxUrlToPath" -> JsObject(fileUriToPath),
        "folderDxUrlToPath" -> JsObject(folderUriToPath)
    )
    FileUtils.writeFileContent(jobMeta.workerPaths.getTaskEnvFile(), json.prettyPrint)
  }

  private def readEnv(): (Map[String, Type],
                          Map[String, (Type, Value)],
                          Map[AddressableFileNode, Path],
                          Map[AddressableFileSource, Path]) = {
    val (schemasJs, inputsJs, filesJs, foldersJs) =
      FileUtils.readFileContent(jobMeta.workerPaths.getTaskEnvFile()).parseJson match {
        case env: JsObject =>
          env
            .getFields("schemas", "localizedInputs", "fileDxUrlToPath", "folderDxUrlToPath") match {
            case Seq(JsObject(schemas),
                     JsObject(inputs),
                     JsObject(filePaths),
                     JsObject(folderPaths)) =>
              (schemas, inputs, filePaths, folderPaths)
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
              val irValue = ValueSerde.deserializeWithType(valueJs, irType)
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
    val folderSourceToPath = foldersJs.map {
      case (uri, JsString(path)) => jobMeta.fileResolver.resolveDirectory(uri) -> Paths.get(path)
      case other                 => throw new Exception(s"unexpected path ${other}")
    }
    (newSchemas, inputs, fileSourceToPath, folderSourceToPath)
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

  /**
    * Collection of files associated with an input value.
    * @param virutalFiles VFiles that have a `contents` attribute - these are
    *                     localized by writing the contents to a file
    * @param realFiles local/remote files
    * @param folders local/remove folders - these are localized by streaming/
    *                downloading the folder recursively
    */
  case class PathValues(virutalFiles: Vector[VFile] = Vector.empty,
                        realFiles: Vector[(VFile, AddressableFileNode)] = Vector.empty,
                        folders: Vector[(VFolder, AddressableFileSource)] = Vector.empty)

  // TODO: it would be nice to extract dx:// links from VString values - this will
  //  happen in the case where the container is a dx file and being passed in as
  //  an input parameter - so that they could be downloaded using dxda. However,
  //  this would also require some way for the downloaded image tarball to be
  //  discovered and loaded. For now, we rely on DockerUtils to download the image
  //  (via DxFileSource, which uses the dx API to download the file).
  private def extractPaths(v: Value): PathValues = {
    def extractArray(values: Vector[Value]): PathValues = {
      values.map(extractPaths).foldLeft(PathValues()) {
        case (pathsAccu, paths) =>
          PathValues(pathsAccu.virutalFiles ++ paths.virutalFiles,
                     pathsAccu.realFiles ++ paths.realFiles,
                     pathsAccu.folders ++ paths.folders)
      }
    }
    v match {
      case file: VFile if file.contents.nonEmpty =>
        PathValues(virutalFiles = Vector(file))
      case file: VFile =>
        val secondaryPaths = extractArray(file.secondaryFiles)
        val fileSource = fileResolver.resolve(file.uri)
        secondaryPaths.copy(realFiles = secondaryPaths.realFiles :+ (file, fileSource))
      case folder: VFolder =>
        val folderSource = fileResolver.resolveDirectory(folder.uri)
        PathValues(folders = Vector((folder, folderSource)))
      case listing: VListing => extractArray(listing.items)
      case VArray(items)     => extractArray(items)
      case VHash(m)          => extractArray(m.values.toVector)
      case _                 => PathValues()
    }
  }

  /**
    * For all input files/directories, determines the local path, whether they need
    * to be streamed or downloaded, and generates the dxda and dxfuse manifests.
    * @return
    * TODO: ensure secondary files are in the same dir as the parent file
    */
  def localizeInputFiles: (Map[String, (Type, Value)],
                           Map[AddressableFileNode, Path],
                           Map[AddressableFileSource, Path],
                           Option[DxdaManifest],
                           Option[DxfuseManifest]) = {

    assert(jobMeta.workerPaths.getInputFilesDir() != jobMeta.workerPaths.getDxfuseMountDir())

    val inputs = getInputsWithDefaults

    case class PathsToLocalize(virtualFiles: Vector[VFile] = Vector.empty,
                               localFiles: Set[(VFile, LocalFileSource)] = Set.empty,
                               filesToStream: Set[(VFile, AddressableFileNode)] = Set.empty,
                               filesToDownload: Set[(VFile, AddressableFileNode)] = Set.empty,
                               localFolders: Set[(VFolder, LocalFileSource)] = Set.empty,
                               foldersToStream: Set[(VFolder, AddressableFileSource)] = Set.empty,
                               foldersToDownload: Set[(VFolder, AddressableFileSource)] = Set.empty)

    def splitFiles(
        files: Vector[(VFile, AddressableFileNode)],
        stream: Boolean
    ): (Set[(VFile, LocalFileSource)],
        Set[(VFile, AddressableFileNode)],
        Set[(VFile, AddressableFileNode)]) = {
      val (local, remote) = files.foldLeft(
          (Set.empty[(VFile, LocalFileSource)], Set.empty[(VFile, AddressableFileNode)])
      ) {
        case (_, (_, fs: LocalFileSource)) if !fs.exists =>
          throw new Exception(s"Local File-type input does not exist: ${fs}")
        case (_, (_, fs: LocalFileSource)) if fs.isDirectory =>
          throw new Exception(s"Local File-type input is a directory: ${fs}")
        case ((local, remote), (file, fs: LocalFileSource)) =>
          // The file is already on the local disk, there is no need to download it.
          (local + (file -> fs), remote)
        case ((local, remote), other) =>
          (local, remote + other)
      }
      if (stream) {
        (local, remote, Set.empty)
      } else {
        (local, Set.empty, remote)
      }
    }

    // splits a Vector of AddressableFileSource associated with an input parameter into
    // (local_files, files_to_stream, files_to_download)
    def splitFolders(
        folders: Vector[(VFolder, AddressableFileSource)],
        stream: Boolean
    ): (Set[(VFolder, LocalFileSource)],
        Set[(VFolder, AddressableFileSource)],
        Set[(VFolder, AddressableFileSource)]) = {
      val (local, remote) = folders.foldLeft(
          (Set.empty[(VFolder, LocalFileSource)], Set.empty[(VFolder, AddressableFileSource)])
      ) {
        case (_, (_, fs: LocalFileSource)) if !fs.exists =>
          throw new Exception(s"Local Directory-type input does not exist: ${fs}")
        case (_, (_, fs: LocalFileSource)) if !fs.isDirectory =>
          throw new Exception(s"Local Directory-type input is a file: ${fs}")
        case ((local, remote), (folder, fs: LocalFileSource)) =>
          // The file is already on the local disk, there is no need to download it.
          (local + (folder -> fs), remote)
        case ((local, remote), other) =>
          (local, remote + other)
      }
      if (stream) {
        (local, remote, Set.empty)
      } else {
        (local, Set.empty, remote)
      }
    }

    val paths = inputs.foldLeft(PathsToLocalize()) {
      case (paths, (name, (_, irValue))) =>
        // extract all the files/directories nested within the input value
        val fileSources = extractPaths(irValue)
        // whether to stream all the files associated with this input
        val stream = streamFiles == StreamFiles.All ||
          (streamFiles == StreamFiles.PerFile && streamFileForInput(name))
        val (localFilesToPath, filesToStream, filesToDownload) =
          splitFiles(fileSources.realFiles, stream)
        val (localFolders, foldersToStream, foldersToDownload) =
          splitFolders(fileSources.folders, stream)
        PathsToLocalize(
            fileSources.virutalFiles,
            paths.localFiles ++ localFilesToPath,
            paths.filesToStream ++ filesToStream,
            paths.filesToDownload ++ filesToDownload,
            paths.localFolders ++ localFolders,
            paths.foldersToStream ++ foldersToStream,
            paths.foldersToDownload ++ foldersToDownload
        )
    }

    // This object handles mapping file sources to local paths and
    // deals with file name collisions in the manner specified by
    // the WDL spec.
    val localizer = SafeLocalizationDisambiguator(
        jobMeta.workerPaths.getInputFilesDir(),
        separateDirsBySource = true,
        createDirs = true,
        disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
        logger = logger
    )

    // localize any virtual files
    logger.traceLimited(s"virtual files = ${paths.virtualFiles}")
    val virtualUriToFile = paths.virtualFiles.map { f =>
      val localPath = localizer.getLocalPath(f.uri, "<virtual>")
      FileUtils.writeFileContent(localPath, f.contents.get)
      f.uri -> f.copy(uri = localPath.toString)
    }.toMap

    // create links for any local files/folders with basename set
    val localFileToPath: Map[AddressableFileNode, Path] = paths.localFiles.map {
      case (file, fs) if file.basename.isDefined =>
        val renamedPath = fs.canonicalPath.getParent.resolve(file.basename.get)
        val localizedPath = localizer.getLocalPath(jobMeta.fileResolver.fromPath(renamedPath))
        Files.createSymbolicLink(localizedPath, fs.canonicalPath)
        fs -> localizedPath
      case (_, fs) => fs -> fs.canonicalPath
    }.toMap
    val localFolderToPath: Map[AddressableFileSource, Path] = paths.localFolders.map {
      case (folder, fs) if folder.basename.isDefined =>
        val renamedPath = fs.canonicalPath.getParent.resolve(folder.basename.get)
        val localizedPath = localizer.getLocalPath(jobMeta.fileResolver.fromPath(renamedPath))
        Files.createSymbolicLink(localizedPath, fs.canonicalPath)
        fs -> localizedPath
      case (_, fs) => fs -> fs.canonicalPath
    }.toMap

    // Build dxda manifest to localize all non-streaming remote files and folders
    logger.traceLimited(s"downloading files = ${paths.filesToDownload}")
    val downloadFileSourceToPath: Map[AddressableFileNode, Path] =
      localizer.getLocalPaths(paths.filesToDownload.toMap.values)
    val downloadFileToPath = downloadFileSourceToPath.collect {
      case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
    }
    val downloadFolderSourceToPath: Map[AddressableFileSource, Path] =
      localizer.getLocalPaths(paths.foldersToDownload.toMap.values)
    val downloadFolderToPath = downloadFolderSourceToPath.collect {
      case (dxFs: DxFolderSource, localPath) => (dxFs.dxProject.id, dxFs.folder) -> localPath
    }
    val dxdaManifest =
      DxdaManifestBuilder(dxApi, logger).apply(downloadFileToPath, downloadFolderToPath)

    // build dxfuse manifest to localize all straming remote files and folders
    logger.traceLimited(s"streaming files = ${paths.filesToStream ++ paths.filesToDownload}")
    val streamingLocalizer = SafeLocalizationDisambiguator(
        jobMeta.workerPaths.getDxfuseMountDir(),
        existingPaths = localizer.getLocalizedPaths,
        separateDirsBySource = true,
        createDirs = false,
        disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
        logger = logger
    )
    val streamFileSourceToPath: Map[AddressableFileNode, Path] =
      streamingLocalizer.getLocalPaths(paths.filesToStream.toMap.values)
    val streamFilesToPaths = streamFileSourceToPath.collect {
      case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
    }
    val streamFolderSourceToPath: Map[AddressableFileSource, Path] =
      streamingLocalizer.getLocalPaths(paths.foldersToStream.toMap.values)
    val streamFoldersToPath = streamFolderSourceToPath.collect {
      case (dxFs: DxFolderSource, localPath) => (dxFs.dxProject.id, dxFs.folder) -> localPath
    }
    val dxfuseManifest = DxfuseManifestBuilder(dxApi, logger)
      .apply(streamFilesToPaths, streamFoldersToPath, jobMeta.workerPaths)

    val fileSourceToPath = localFileToPath ++ downloadFileSourceToPath ++ streamFileSourceToPath
    val folderSourceToPath = localFolderToPath ++ downloadFolderSourceToPath ++ streamFolderSourceToPath
    val uriToPath: Map[String, String] =
      (fileSourceToPath ++ folderSourceToPath).map {
        case (dxFs: DxFileSource, path)       => dxFs.address -> path.toString
        case (dxFs: DxFolderSource, path)     => dxFs.address -> path.toString
        case (localFs: LocalFileSource, path) => localFs.originalPath.toString -> path.toString
        case (other, _)                       => throw new RuntimeException(s"unsupported file source ${other}")
      }.toMap

    // Replace the URIs with local file paths
    def pathTranslator(v: Value, t: Option[Type], optional: Boolean): Option[PathValue] = {
      def translateUri(uri: String): Option[String] = {
        uriToPath.get(uri) match {
          case Some(localPath)  => Some(localPath)
          case None if optional => None
          case _                => throw new Exception(s"Did not localize file ${uri}")
        }
      }
      (t, v) match {
        case (_, f: VFile) if f.contents.nonEmpty =>
          Some(virtualUriToFile(f.uri))
        case (_, f: VFile) =>
          translateUri(f.uri).map { newUri =>
            val newSecondaryFiles =
              f.secondaryFiles.flatMap(pathTranslator(_, None, optional = false))
            f.copy(uri = newUri, secondaryFiles = newSecondaryFiles)
          }
        case (Some(TFile), VString(uri)) => translateUri(uri).map(VFile(_))
        case (_, f: VFolder)             => translateUri(f.uri).map(newUri => f.copy(uri = newUri))
        case (_, l: VListing) =>
          val newListing = l.items.map(pathTranslator(_, None, optional = optional))
          if (newListing.forall(_.isEmpty)) {
            None
          } else if (newListing.exists(_.isEmpty)) {
            val notLocalized = l.items.zip(newListing).collect {
              case (old, None) => old
            }
            throw new Exception(
                s"some listing entries were not localized: ${notLocalized.mkString(",")}"
            )
          } else {
            Some(l.copy(items = newListing.flatten))
          }
        case _ => None
      }
    }

    val localizedInputs = inputs.view.mapValues {
      case (t, v) => (t, Value.transform(v, Some(t), pathTranslator))
    }.toMap

    (localizedInputs, fileSourceToPath, folderSourceToPath, dxdaManifest, dxfuseManifest)
  }

  protected def getSchemas: Map[String, TSchema]

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

    val (localizedInputs, fileSourceToPath, folderSourceToPath, dxdaManifest, dxfuseManifest) =
      localizeInputFiles

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

    writeEnv(getSchemas, localizedInputs, fileSourceToPath, folderSourceToPath)
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
    val (schemas, localizedInputs, fileSourceToPath, folderSourceToPath) =
      readEnv()
    logger.traceLimited(s"InstantiateCommand, env = ${localizedInputs}")
    // evaluate the command block and write the command script
    val updatedInputs = writeCommandScript(localizedInputs)
    // write the updated env to disk
    writeEnv(schemas, updatedInputs, fileSourceToPath, folderSourceToPath)
  }

  /**
    * Evaluates the outputs of the task, uploads and de-localizes output files,
    * and write outputs to the meta file.
    * @param localizedInputs the job inputs, with files localized to the worker
    */
  protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)]

  private lazy val virtualOutputDir =
    Files.createTempDirectory(jobMeta.workerPaths.getOutputFilesDir(ensureExists = true), "virtual")

  // TODO: handle secondary files - put files in same output folder as the primary
  // TODO: handle basename
  private def extractOutputFiles(
      name: String,
      v: Value,
      t: Type
  ): (Vector[(VFile, AddressableFileNode)], Vector[(VFolder, AddressableFileSource)]) = {
    def getLocalFileSource(varName: String,
                           fs: AddressableFileSource,
                           optional: Boolean): Option[LocalFileSource] = {
      fs match {
        case localFs: LocalFileSource if optional && !Files.exists(localFs.canonicalPath) =>
          // ignore optional, non-existent files
          None
        case localFs: LocalFileSource if !Files.exists(localFs.canonicalPath) =>
          throw new Exception(
              s"required output file ${varName} does not exist at ${localFs.canonicalPath}"
          )
        case localFs: LocalFileSource =>
          Some(localFs)
        case other =>
          throw new RuntimeException(s"${varName} specifies non-local file ${other}")
      }
    }

    def inner(
        innerName: String,
        innerValue: Value,
        innerType: Type,
        optional: Boolean
    ): (Vector[(VFile, LocalFileSource)], Vector[(VFolder, LocalFileSource)]) = {
      (innerType, innerValue) match {
        case (TOptional(_), VNull) =>
          (Vector.empty, Vector.empty)
        case (TOptional(optType), _) =>
          inner(innerName, innerValue, optType, optional = true)
        case (TFile, f: VFile) if f.contents.isDefined =>
          val path = f.basename
            .map(virtualOutputDir.resolve)
            .getOrElse(Files.createTempFile(virtualOutputDir, "virtual", ""))
          val resolved = FileUtils.writeFileContent(path, f.contents.get)
          (Vector((f, fileResolver.fromPath(resolved))), Vector.empty)
        case (TFile, f: VFile) =>
          (getLocalFileSource(innerName, fileResolver.resolve(f.uri), optional)
             .map(f -> _)
             .toVector,
           Vector.empty)
        case (TFile, VString(path)) =>
          (getLocalFileSource(innerName, fileResolver.resolve(path), optional)
             .map(VFile(path) -> _)
             .toVector,
           Vector.empty)
        case (TDirectory, f: VFolder) =>
          (Vector.empty,
           getLocalFileSource(innerName, fileResolver.resolveDirectory(f.uri), optional)
             .map(f -> _)
             .toVector)
        case (TDirectory, f: VListing) =>
          val (files, folders) = f.items.zipWithIndex.map {
            case (p: PathValue, index) =>
              inner(s"${innerName}/${index}", p, TDirectory, optional)
          }.unzip
          (files.flatten, folders.flatten)
        case (TArray(_, nonEmpty), VArray(array)) if nonEmpty && array.isEmpty =>
          throw new Exception(s"Non-empty array ${name} has empty value")
        case (TArray(elementType, _), VArray(array)) =>
          val (files, folders) = array.zipWithIndex.map {
            case (element, index) =>
              inner(s"${innerName}[${index}]", element, elementType, optional = false)
          }.unzip
          (files.flatten, folders.flatten)
        case (TSchema(name, memberTypes), VHash(members)) =>
          val (files, folders) = memberTypes
            .map {
              case (key, t) =>
                (t, members.get(key)) match {
                  case (TOptional(_), None) =>
                    (Vector.empty, Vector.empty)
                  case (_, None) =>
                    throw new Exception(s"missing non-optional member ${key} of struct ${name}")
                  case (_, Some(v)) =>
                    inner(s"${name}.${key}", v, t, optional = false)
                }
            }
            .toVector
            .unzip
          (files.flatten, folders.flatten)
        case (THash, VHash(members)) =>
          val (files, folders) = members.toVector.map {
            case (key, value) =>
              val fileSources = extractPaths(value)
              val files = fileSources.realFiles.flatMap {
                case (file, fs) =>
                  getLocalFileSource(s"${innerName}.${key}", fs, optional = true)
                    .map(file -> _)
                    .toVector
              }
              val folders = fileSources.folders.flatMap {
                case (folder, fs) =>
                  getLocalFileSource(s"${innerName}.${key}", fs, optional = true)
                    .map(folder -> _)
                    .toVector
              }
              (files, folders)
          }.unzip
          (files.flatten, folders.flatten)
        case _ =>
          (Vector.empty, Vector.empty)
      }
    }

    inner(name, v, t, optional = false)
  }

  def epilog(): Unit = {
    if (logger.isVerbose) {
      trace(s"Epilog debugLevel=${logger.traceLevel}")
      printDirTree()
    }
    val (_, localizedInputs, fileSourceToPath, _) = readEnv()
    val localizedOutputs = evaluateOutputs(localizedInputs)

    // extract files from the outputs
    val (localOutputFiles, localOutputFolders) = localizedOutputs.map {
      case (name, (irType, irValue)) => extractOutputFiles(name, irValue, irType)
    }.unzip

    // Build a map of all the string values in the output values that might
    // map to the same (absolute) local path. Some of the outputs may be files
    // that were inputs (in `fileSourceToPath`) - these do not need to be
    // re-uploaded. The `localPath`s will be the same but the `originalPath`s
    // may be different.
    val inputAddresses: Set[String] = fileSourceToPath.keySet.map(_.address)
    val inputPaths = fileSourceToPath.values.toSet
    val delocalizingFileUriToPath: Map[VFile, Path] = localOutputFiles.flatten.collect {
      case (file, local: LocalFileSource)
          if !(
              inputAddresses.contains(file.uri) || inputPaths.contains(local.canonicalPath)
          ) =>
        file -> local.canonicalPath
    }.toMap
    val delocalizingFolderUriToPath: Map[VFolder, Path] = localOutputFolders.flatten.collect {
      case (folder, local: LocalFileSource)
          if !(
              inputAddresses.contains(folder.uri) || inputPaths.contains(local.canonicalPath)
          ) =>
        folder -> local.canonicalPath
    }.toMap

    // upload the files, and map their local paths to their remote URIs
    val delocalizedPathToUri: Map[Path, String] = {
      val (dxFiles, dxFolders) = if (jobMeta.useManifests) {
        // if using manifests, we need to upload the files directly to the project
        (fileUploader.uploadFilesWithDestination(
             delocalizingFileUriToPath.toSet.map {
               case (file, path) =>
                 val basename = file.basename.getOrElse(path.getFileName.toString)
                 path -> s"${jobMeta.manifestFolder}/${basename}"
             }.toMap,
             wait = waitOnUpload
         ),
         fileUploader.uploadDirectoriesWithDestination(
             delocalizingFolderUriToPath.toSet.map {
               case (folder, path) =>
                 val basename = folder.basename.getOrElse(path.getFileName.toString)
                 path -> s"${jobMeta.manifestFolder}/${basename}"
             }.toMap,
             wait = waitOnUpload
         ))
      } else {
        (fileUploader.uploadFiles(files.toSet, wait = waitOnUpload),
         fileUploader.uploadDirectories(dirs.toSet, wait = waitOnUpload))
      }
      val fileUris = dxFiles.map {
        case (path, dxFile) => path -> dxFile.asUri
      }
      val folderUris = dxFolders.map {
        case (path, (projectId, folder)) => path -> s"${DxPath.DxUriPrefix}${projectId}:${folder}"
      }
      fileUris ++ folderUris
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
        .orElse(delocalizingUriToPath.get(value) match {
          case Some(path) => delocalizedPathToUri.get(path)
          case _          => None
        })
    }

    // Replace the URIs with remote file paths
    def pathTranslator(v: Value, t: Option[Type], optional: Boolean): Option[PathValue] = {
      def translateUri(uri: String): Option[String] = {
        resolveFileValue(uri) match {
          case Some(localPath)  => Some(localPath)
          case None if optional => None
          case _                => throw new Exception(s"Did not localize file ${uri}")
        }
      }
      (t, v) match {
        case (_, f: VFile) =>
          translateUri(f.uri).map { newUri =>
            val newSecondaryFiles =
              f.secondaryFiles.flatMap(pathTranslator(_, None, optional = false))
            f.copy(uri = newUri, secondaryFiles = newSecondaryFiles)
          }
        case (Some(TFile), VString(uri)) => translateUri(uri).map(VFile(_))
        case (_, f: VFolder)             => translateUri(f.uri).map(newUri => f.copy(uri = newUri))
        case (_, l: VListing) =>
          val newListing = l.items.map(pathTranslator(_, None, optional = optional))
          if (newListing.forall(_.isEmpty)) {
            None
          } else if (newListing.exists(_.isEmpty)) {
            val notLocalized = l.items.zip(newListing).collect {
              case (old, None) => old
            }
            throw new Exception(
                s"some listing entries were not localized: ${notLocalized.mkString(",")}"
            )
          } else {
            Some(l.copy(items = newListing.flatten))
          }
        case _ => None
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
