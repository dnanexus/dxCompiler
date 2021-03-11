package dx.executor

import java.nio.file.{Files, Path, Paths}
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
import dx.util.protocols.{DxArchiveFolderSource, DxFileSource, DxFolderSource}
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
                            streamFiles: StreamFiles.StreamFiles,
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
                       archiveSourceToPath: Map[AddressableFileNode, Path],
                       folderSourceToPath: Map[AddressableFileSource, Path]): Unit = {
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
    val fileUriToPath: Map[String, JsValue] = fileSourceToPath.map {
      case (fileSource: AddressableFileSource, path) =>
        fileSource.address -> JsString(path.toString)
      case (other, _) =>
        throw new RuntimeException(s"Can only serialize an AddressableFileSource, not ${other}")
    }
    val archiveUriToPath: Map[String, JsValue] = archiveSourceToPath.map {
      case (archiveSource: AddressableFileSource, path) =>
        archiveSource.address -> JsString(path.toString)
    }
    val folderUriToPath: Map[String, JsValue] = folderSourceToPath.map {
      case (folderSource: AddressableFileSource, path) =>
        folderSource.address -> JsString(path.toString)
    }
    val json = JsObject(
        "schemas" -> JsObject(newSchemasJs),
        "localizedInputs" -> JsObject(inputsJs),
        "fileDxUrlToPath" -> JsObject(fileUriToPath),
        "archiveDxUrlToPath" -> JsObject(archiveUriToPath),
        "folderDxUrlToPath" -> JsObject(folderUriToPath)
    )
    FileUtils.writeFileContent(jobMeta.workerPaths.getTaskEnvFile(), json.prettyPrint)
  }

  private def readEnv(): (Map[String, Type],
                          Map[String, (Type, Value)],
                          Map[AddressableFileNode, Path],
                          Map[AddressableFileNode, Path],
                          Map[AddressableFileSource, Path]) = {
    val (schemasJs, inputsJs, filesJs, archivesJs, foldersJs) =
      FileUtils.readFileContent(jobMeta.workerPaths.getTaskEnvFile()).parseJson match {
        case env: JsObject =>
          env.getFields("schemas",
                        "localizedInputs",
                        "fileDxUrlToPath",
                        "archiveDxUrlToPath",
                        "folderDxUrlToPath") match {
            case Seq(JsObject(schemas),
                     JsObject(inputs),
                     JsObject(filePaths),
                     JsObject(archivePaths),
                     JsObject(folderPaths)) =>
              (schemas, inputs, filePaths, archivePaths, folderPaths)
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
    val archiveSourceToPath = archivesJs.map {
      case (uri, JsString(path)) => jobMeta.fileResolver.resolve(uri) -> Paths.get(path)
      case other                 => throw new Exception(s"unexpected path ${other}")
    }
    val folderSourceToPath = foldersJs.map {
      case (uri, JsString(path)) => jobMeta.fileResolver.resolveDirectory(uri) -> Paths.get(path)
      case other                 => throw new Exception(s"unexpected path ${other}")
    }
    (newSchemas, inputs, fileSourceToPath, archiveSourceToPath, folderSourceToPath)
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

  // TODO: it would be nice to extract dx:// links from VString values - this will
  //  happen in the case where the container is a dx file and being passed in as
  //  an input parameter - so that they could be downloaded using dxda. However,
  //  this would also require some way for the downloaded image tarball to be
  //  discovered and loaded. For now, we rely on DockerUtils to download the image
  //  (via DxFileSource, which uses the dx API to download the file).
  private def extractPaths(
      v: Value
  ): (Vector[AddressableFileNode], Vector[AddressableFileNode], Vector[AddressableFileSource]) = {
    def extractArray(
        values: Vector[Value]
    ): (Vector[AddressableFileNode], Vector[AddressableFileNode], Vector[AddressableFileSource]) = {
      val (files, archives, folders) = values.map(extractPaths).unzip3
      (files.flatten, archives.flatten, folders.flatten)
    }
    v match {
      case f: VFile =>
        val primaryFile = fileResolver.resolve(f.uri)
        val (secFiles, secArchives, secFolders) = extractArray(f.secondaryFiles)
        (secFiles :+ primaryFile, secArchives, secFolders)
      case f: VFolder =>
        (Vector.empty, Vector.empty, Vector(fileResolver.resolveDirectory(f.uri)))
      case a: VArchive =>
        (Vector.empty, Vector(fileResolver.resolve(a.uri)), Vector.empty)
      case l: VListing   => extractArray(l.listing)
      case VArray(items) => extractArray(items)
      case VHash(m)      => extractArray(m.values.toVector)
      case _             => (Vector.empty, Vector.empty, Vector.empty)
    }
  }

  /**
    * For all input files/directories, determines the local path, whether they need
    * to be streamed or downloaded, and generates the dxda and dxfuse manifests.
    * @return
    * TODO: handle file/archive/folder with basename
    * TODO: handle file with contents
    * TODO: basedir for paths in listing
    * TODO: ensure secondary files are in the same dir as the parent file
    */
  def localizeInputFiles: (Map[String, (Type, Value)],
                           Map[AddressableFileNode, Path],
                           Map[AddressableFileNode, Path],
                           Map[AddressableFileSource, Path],
                           Option[DxdaManifest],
                           Option[DxfuseManifest]) = {

    assert(jobMeta.workerPaths.getInputFilesDir() != jobMeta.workerPaths.getDxfuseMountDir())

    val inputs = getInputsWithDefaults

    case class Paths(localFilesToPath: Map[AddressableFileNode, Path] = Map.empty,
                     filesToStream: Set[AddressableFileNode] = Set.empty,
                     filesToDownload: Set[AddressableFileNode] = Set.empty,
                     localArchivesToPath: Map[AddressableFileNode, Path] = Map.empty,
                     archivesToStream: Set[AddressableFileNode] = Set.empty,
                     archivesToDownload: Set[AddressableFileNode] = Set.empty,
                     localFoldersToPath: Map[AddressableFileSource, Path] = Map.empty,
                     foldersToStream: Set[AddressableFileSource] = Set.empty) {
      lazy val localPaths: Set[Path] =
        (localFilesToPath.values ++ localArchivesToPath.values ++ localFoldersToPath.values).toSet
    }

    def splitFiles(
        name: String,
        files: Vector[AddressableFileNode]
    ): (Map[AddressableFileNode, Path], Set[AddressableFileNode], Set[AddressableFileNode]) = {
      val (local, remote) = files.foldLeft(
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
        (local, remote, Set.empty)
      } else {
        (local, Set.empty, remote)
      }
    }

    val paths = inputs.foldLeft(Paths()) {
      case (paths, (name, (_, irValue))) =>
        val (files, archives, folders) = extractPaths(irValue)
        val (localFilesToPath, filesToStream, filesToDownload) = splitFiles(name, files)
        val (localArchivesToPath, archivesToStream, archivesToDownload) = splitFiles(name, archives)
        val (localFolders, remoteFolders) = folders.foldLeft(Map.empty[AddressableFileSource, Path],
                                                             Set.empty[AddressableFileSource]) {
          case ((local, remote), fs: LocalFileSource) if fs.isDirectory =>
            (local + (fs -> fs.canonicalPath), remote)
          case ((local, remote), other) =>
            (local, remote + other)
        }
        Paths(
            paths.localFilesToPath ++ localFilesToPath,
            paths.filesToStream ++ filesToStream,
            paths.filesToDownload ++ filesToDownload,
            paths.localArchivesToPath ++ localArchivesToPath,
            paths.archivesToStream ++ archivesToStream,
            paths.archivesToDownload ++ archivesToDownload,
            paths.localFoldersToPath ++ localFolders,
            paths.foldersToStream ++ remoteFolders
        )
    }

    // build dxda and/or dxfuse manifests
    // We use a SafeLocalizationDisambiguator to determine the local path and deal
    // with file name collisions in the manner specified by the WDL spec. We set
    // separateDirsBySource = true because, when creating archives, we append one
    // directory at a time to the archive and then delete the source files, so the
    // smaller each append operation is, the less disk overhead is required.

    logger.traceLimited(s"downloading files = ${paths.filesToDownload ++ paths.archivesToDownload}")
    val downloadLocalizer =
      SafeLocalizationDisambiguator(
          jobMeta.workerPaths.getInputFilesDir(),
          existingPaths = paths.localPaths,
          separateDirsBySource = true,
          disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs
      )
    val downloadFileSourceToPath =
      paths.filesToDownload.map(fs => fs -> downloadLocalizer.getLocalPath(fs)).toMap
    val downloadArchiveSourceToPath =
      paths.archivesToDownload.map(fs => fs -> downloadLocalizer.getLocalPath(fs)).toMap
    // write the manifest for dxda, if there are files to download
    val dxdaManifest = DxdaManifestBuilder(dxApi)
      .apply((downloadFileSourceToPath ++ downloadArchiveSourceToPath).collect {
        case (dxFs: DxFileSource, localPath)          => dxFs.dxFile -> localPath
        case (dxFs: DxArchiveFolderSource, localPath) => dxFs.dxFileSource.dxFile -> localPath
      })

    logger.traceLimited(s"streaming files = ${paths.filesToStream ++ paths.filesToDownload}")
    val streamingLocalizer =
      SafeLocalizationDisambiguator(
          jobMeta.workerPaths.getDxfuseMountDir(),
          existingPaths = paths.localPaths,
          separateDirsBySource = true,
          disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs
      )
    val streamFileSourceToPath =
      paths.filesToStream.map(fs => fs -> streamingLocalizer.getLocalPath(fs)).toMap
    val streamArchiveSourceToPath =
      paths.archivesToStream.map(fs => fs -> streamingLocalizer.getLocalPath(fs)).toMap
    val dxFilesToPaths = (streamFileSourceToPath ++ streamArchiveSourceToPath).collect {
      case (dxFs: DxFileSource, localPath)          => dxFs.dxFile -> localPath
      case (dxFs: DxArchiveFolderSource, localPath) => dxFs.dxFileSource.dxFile -> localPath
    }
    val streamFolderSourceToPath: Map[AddressableFileSource, Path] =
      paths.foldersToStream.map(fs => fs -> streamingLocalizer.getLocalPath(fs)).toMap
    val dxFoldersToPath = streamFolderSourceToPath.collect {
      case (dxFs: DxFolderSource, localPath) => (dxFs.dxProject.id, dxFs.folder) -> localPath
    }
    // write the manifest for dxfuse, if there are files to stream
    val dxfuseManifest = DxfuseManifestBuilder(dxApi)
      .apply(dxFilesToPaths, dxFoldersToPath, jobMeta.workerPaths)

    val fileSourceToPath = paths.localFilesToPath ++ downloadFileSourceToPath ++ streamFileSourceToPath
    val archiveSourceToPath = paths.localArchivesToPath ++ downloadArchiveSourceToPath ++ streamArchiveSourceToPath
    val folderSourceToPath = paths.localFoldersToPath ++ streamFileSourceToPath

    val uriToPath: Map[String, String] =
      (fileSourceToPath ++ archiveSourceToPath ++ folderSourceToPath).map {
        case (dxFs: DxFileSource, path)          => dxFs.address -> path.toString
        case (dxFs: DxFolderSource, path)        => dxFs.address -> path.toString
        case (dxFs: DxArchiveFolderSource, path) => dxFs.address -> path.toString
        case (localFs: LocalFileSource, p2)      => localFs.originalPath.toString -> p2.toString
        case other                               => throw new RuntimeException(s"unsupported file source ${other}")
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
        case (_, f: VFile) =>
          translateUri(f.uri).map { newUri =>
            val newSecondaryFiles =
              f.secondaryFiles.flatMap(pathTranslator(_, None, optional = false))
            f.copy(uri = newUri, secondaryFiles = newSecondaryFiles)
          }
        case (Some(TFile), VString(uri)) => translateUri(uri).map(VFile(_))
        case (_, f: VFolder)             => translateUri(f.uri).map(newUri => f.copy(uri = newUri))
        case (_, a: VArchive)            => translateUri(a.uri).map(newUri => a.copy(uri = newUri))
        case (_, l: VListing) =>
          val newListing = l.listing.map(pathTranslator(_, None, optional = optional))
          if (newListing.forall(_.isEmpty)) {
            None
          } else if (newListing.exists(_.isEmpty)) {
            val notLocalized = l.listing.zip(newListing).collect {
              case (old, None) => old
            }
            throw new Exception(
                s"some listing entries were not localized: ${notLocalized.mkString(",")}"
            )
          } else {
            Some(l.copy(listing = newListing.flatten))
          }
        case _ => None
      }
    }

    val localizedInputs = inputs.view.mapValues {
      case (t, v) => (t, Value.transform(v, Some(t), pathTranslator))
    }.toMap

    (localizedInputs,
     fileSourceToPath,
     archiveSourceToPath,
     folderSourceToPath,
     dxdaManifest,
     dxfuseManifest)
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

    val (localizedInputs,
         fileSourceToPath,
         archiveSourceToPath,
         folderSourceToPath,
         dxdaManifest,
         dxfuseManifest) = localizeInputFiles

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

    writeEnv(getSchemas, localizedInputs, fileSourceToPath, archiveSourceToPath, folderSourceToPath)
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
    val (schemas, localizedInputs, fileSourceToPath, archiveSourceToPath, folderSourceToPath) =
      readEnv()
    logger.traceLimited(s"InstantiateCommand, env = ${localizedInputs}")
    // unpack any archives and update the paths
    val archiveSourceToUnpackedPath = archiveSourceToPath.map {
      case (fs, path) =>
        // TODO: use the basename here - for now just using the filename without the extension
        val folder = FileUtils.changeFirstFileExt(fs.name, Vector(".tgz", ".tar.gz", ".tar"))
        val tarOpts = if (fs.name.endsWith("gz")) "-xz" else "-x"
        try {
          SysUtils.execCommand(
              s"mkdir ${folder} && tar ${tarOpts} -f ${path} -C ${folder}"
          )
        } finally {
          path.toFile.delete()
        }
        fs -> path.getParent.resolve(folder)
    }
    // evaluate the command block and write the command script
    val updatedInputs = writeCommandScript(localizedInputs)
    // write the updated env to disk
    writeEnv(schemas,
             updatedInputs,
             fileSourceToPath,
             archiveSourceToUnpackedPath,
             folderSourceToPath)
  }

  /**
    * Evaluates the outputs of the task, uploads and de-localizes output files,
    * and write outputs to the meta file.
    * @param localizedInputs the job inputs, with files localized to the worker
    */
  protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)]

  // TODO: handle secondary files - put files in same output folder as the primary
  // TODO: listing basedir
  // TODO: handle basename
  // TODO: handle file with contents
  private def extractOutputFiles(name: String, v: Value, t: Type): Vector[AddressableFileSource] = {
    def getFileNode(varName: String,
                    fs: AddressableFileSource,
                    optional: Boolean): Vector[AddressableFileSource] = {
      fs match {
        case localFs: LocalFileSource if optional && !Files.exists(localFs.canonicalPath) =>
          // ignore optional, non-existent files
          Vector.empty
        case localFs: LocalFileSource if !Files.exists(localFs.canonicalPath) =>
          throw new Exception(
              s"required output file ${varName} does not exist at ${localFs.canonicalPath}"
          )
        case localFs: LocalFileSource =>
          Vector(localFs)
        case other =>
          throw new RuntimeException(s"${varName} specifies non-local file ${other}")
      }
    }

    def inner(innerName: String,
              innerValue: Value,
              innerType: Type,
              optional: Boolean): Vector[AddressableFileSource] = {
      (innerType, innerValue) match {
        case (TOptional(_), VNull) =>
          Vector.empty
        case (TOptional(optType), _) =>
          inner(innerName, innerValue, optType, optional = true)
        case (TFile, f: VFile) =>
          getFileNode(innerName, fileResolver.resolve(f.uri), optional)
        case (TFile, VString(path)) =>
          getFileNode(innerName, fileResolver.resolve(path), optional)
        case (TDirectory, f: VFolder) =>
          getFileNode(innerName, fileResolver.resolveDirectory(f.uri), optional)
        case (TDirectory, f: VListing) =>
          f.listing.zipWithIndex.flatMap {
            case (p: PathValue, index) =>
              inner(s"${innerName}/${index}", p, TDirectory, optional)
          }
        case (TArray(_, nonEmpty), VArray(array)) if nonEmpty && array.isEmpty =>
          throw new Exception(s"Non-empty array ${name} has empty value")
        case (TArray(elementType, _), VArray(array)) =>
          array.zipWithIndex.flatMap {
            case (element, index) =>
              inner(s"${innerName}[${index}]", element, elementType, optional = false)
          }
        case (TSchema(name, memberTypes), VHash(members)) =>
          memberTypes.flatMap {
            case (key, t) =>
              (t, members.get(key)) match {
                case (TOptional(_), None) =>
                  Vector.empty
                case (_, None) =>
                  throw new Exception(s"missing non-optional member ${key} of struct ${name}")
                case (_, Some(v)) =>
                  inner(s"${name}.${key}", v, t, optional = false)
              }
          }.toVector
        case (THash, VHash(members)) =>
          members.flatMap {
            case (key, value) =>
              val (files, _, folders) = extractPaths(value)
              (files ++ folders).flatMap(fs =>
                getFileNode(s"${innerName}.${key}", fs, optional = true)
              )
          }.toVector
        case _ =>
          Vector.empty
      }
    }

    inner(name, v, t, optional = false)
  }

  def epilog(): Unit = {
    if (logger.isVerbose) {
      trace(s"Epilog debugLevel=${logger.traceLevel}")
      printDirTree()
    }
    val (_, localizedInputs, fileSourceToPath, _, _) = readEnv()
    val localizedOutputs = evaluateOutputs(localizedInputs)

    // extract files from the outputs
    val localOutputFileSources = localizedOutputs.flatMap {
      case (name, (irType, irValue)) => extractOutputFiles(name, irValue, irType)
    }.toVector

    // Build a map of all the string values in the output values that might
    // map to the same (absolute) local path. Some of the outputs may be files
    // that were inputs (in `fileSourceToPath`) - these do not need to be
    // re-uploaded. The `localPath`s will be the same but the `originalPath`s
    // may be different.
    val inputAddresses: Set[String] = fileSourceToPath.keySet.map(_.address)
    val inputPaths = fileSourceToPath.values.toSet
    val delocalizingValueToPath: Map[String, Path] = localOutputFileSources.collect {
      case local: LocalFileSource
          if !(
              inputAddresses.contains(local.address) || inputPaths.contains(local.canonicalPath)
          ) =>
        local.address -> local.canonicalPath
    }.toMap

    // TODO: also do this for archives and folders

    // upload the files, and map their local paths to their remote URIs
    val delocalizedPathToUri: Map[Path, String] = {
      val dxFiles = if (jobMeta.useManifests) {
        // if using manifests, we need to upload the files directly to the project
        fileUploader.upload(delocalizingValueToPath.values.map { path =>
          path -> s"${jobMeta.manifestFolder}/${path.getFileName.toString}"
        }.toMap)
      } else {
        fileUploader.upload(delocalizingValueToPath.values.toSet)
      }
      dxFiles.map {
        case (path, dxFile) => path -> dxFile.asUri
      }
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
        case (_, a: VArchive)            => translateUri(a.uri).map(newUri => a.copy(uri = newUri))
        case (_, l: VListing) =>
          val newListing = l.listing.map(pathTranslator(_, None, optional = optional))
          if (newListing.forall(_.isEmpty)) {
            None
          } else if (newListing.exists(_.isEmpty)) {
            val notLocalized = l.listing.zip(newListing).collect {
              case (old, None) => old
            }
            throw new Exception(
                s"some listing entries were not localized: ${notLocalized.mkString(",")}"
            )
          } else {
            Some(l.copy(listing = newListing.flatten))
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
