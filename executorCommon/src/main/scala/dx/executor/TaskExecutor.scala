package dx.executor

import java.nio.file.{Files, Path, Paths}
import dx.api.{DxFile, DxJob, FileUpload, InstanceTypeRequest}
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
import dx.util.{
  AddressableFileNode,
  AddressableFileSource,
  Enum,
  FileSource,
  FileUtils,
  LocalFileSource,
  SafeLocalizationDisambiguator,
  StringFileNode,
  SysUtils,
  TraceLevel
}
import dx.util.protocols.DxFileSource
import spray.json._

object TaskAction extends Enum {
  type TaskAction = Value
  val CheckInstanceType, LocalizeInputs, FinalizeInputs, InstantiateCommand, DelocalizeOutputs,
      Relaunch = Value
}

object TaskExecutor {
  val MaxDisambiguationDirs: Int = 5000
  val MaxConcurrentUploads: Int = 8
  val InitialInputsFile = "serialized_initial_inputs.json"
  val LocalizedInputsFile = "serialized_localized_inputs.json"
  val LocalizedInputsWithCommandBindingsFile =
    "serialized_localized_inputs_with_command_bindings.json"
  val UriToPathFile = "serialized_uri_to_path.json"
  val PathToUriFile = "serialized_path_to_uri.json"
  val LocalizerFile = "serialized_localizer.json"
}

abstract class TaskExecutor(jobMeta: JobMeta,
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
    logger.trace("Computing instance type to request")
    val instanceTypeRequest: InstanceTypeRequest = getInstanceTypeRequest
    logger.traceLimited(s"Requesting instance type: $instanceTypeRequest")
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
    trace(s"Requested instance type: ${requestedInstanceType}")
    val currentInstanceType = jobMeta.instanceType.getOrElse(
        throw new Exception(s"Cannot get instance type for job ${jobMeta.jobId}")
    )
    trace(s"Current instance type: ${currentInstanceType}")
    val isSufficient =
      try {
        jobMeta.instanceTypeDb.matchesOrExceedes(currentInstanceType, requestedInstanceType)
      } catch {
        case ex: Throwable =>
          logger.warning("Error comparing current and requested instance types",
                         exception = Some(ex))
          false
      }
    if (isSufficient) {
      trace(s"The current instance type is sufficient")
    } else {
      trace(
          s"The current instance type is not sufficient; relaunching job with a larger instance type"
      )
    }
    isSufficient
  }

  // These functions write/read the inputs and mappings between file
  // inputs and their local paths on disk. This is used to preserve
  // state between different phases of task execution, so we don't
  // have to re-evaluate expressions every time.

  private def writeState(fileName: String, jsValue: JsValue): Unit = {
    val path = jobMeta.workerPaths.getMetaDir(ensureExists = true).resolve(fileName)
    FileUtils.writeFileContent(path, jsValue.prettyPrint)
  }

  private def readState(fileName: String): JsValue = {
    val path = jobMeta.workerPaths.getMetaDir().resolve(fileName)
    FileUtils.readFileContent(path, mustExist = true).parseJson
  }

  private def serializeInputs(fileName: String,
                              schemas: Map[String, Type],
                              inputs: Map[String, (Type, Value)]): Unit = {
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
    val jsValue = JsObject(
        "schemas" -> JsObject(newSchemasJs),
        "inputs" -> JsObject(inputsJs)
    )
    writeState(fileName, jsValue)
  }

  private def deserializeInputs(
      fileName: String
  ): (Map[String, Type], Map[String, (Type, Value)]) = {
    readState(fileName) match {
      case env: JsObject =>
        env.getFields("schemas", "inputs") match {
          case Seq(JsObject(schemaJs), JsObject(inputJs)) =>
            val schemas = TypeSerde.deserializeSchemas(schemaJs)
            inputJs.foldLeft((schemas, Map.empty[String, (Type, Value)])) {
              case ((schemaAccu, paramAccu), (key, obj: JsObject)) =>
                obj.getFields("type", "value") match {
                  case Seq(typeJs, valueJs) =>
                    val (irType, newSchemas) = TypeSerde.deserialize(typeJs, schemaAccu)
                    val irValue = ValueSerde.deserializeWithType(valueJs, irType, key)
                    (newSchemas, paramAccu + (key -> (irType, irValue)))
                  case _ =>
                    throw new Exception(s"Malformed environment serialized to disk ${fileName}")
                }
              case _ =>
                throw new Exception(s"Malformed environment serialized to disk ${fileName}")
            }
          case _ =>
            throw new Exception(s"Malformed environment serialized to disk ${fileName}")
        }
      case _ =>
        throw new Exception(s"Malformed environment serialized to disk ${fileName}")
    }
  }

  private def serializeUriToPath(uriToPath: Map[String, Path]): Unit = {
    val uriToPathJs: Map[String, JsValue] = uriToPath.map {
      case (uri, path) => uri -> JsString(path.toString)
    }
    writeState(TaskExecutor.UriToPathFile, JsObject(uriToPathJs))
  }

  private def deserializeUriToPath(): Map[String, Path] = {
    readState(TaskExecutor.UriToPathFile) match {
      case JsObject(uriToPath) =>
        uriToPath.map {
          case (uri, JsString(path)) => uri -> Paths.get(path)
          case other                 => throw new Exception(s"unexpected path ${other}")
        }
      case _ =>
        throw new Exception(
            s"Malformed environment serialized to disk ${TaskExecutor.UriToPathFile}"
        )
    }
  }

  private def serializePathToUri(pathToUri: Map[Path, String]): Unit = {
    val pathToUriJs: Map[String, JsValue] = pathToUri.map {
      case (path, uri) => path.toString -> JsString(uri)
    }
    writeState(TaskExecutor.PathToUriFile, JsObject(pathToUriJs))
  }

  private def deserializePathToUri(): Map[Path, String] = {
    readState(TaskExecutor.PathToUriFile) match {
      case JsObject(pathToUri) =>
        pathToUri.map {
          case (path, JsString(uri)) => Paths.get(path) -> uri
          case other                 => throw new Exception(s"unexpected uri ${other}")
        }
      case _ =>
        throw new Exception(
            s"Malformed environment serialized to disk ${TaskExecutor.PathToUriFile}"
        )
    }
  }

  private def serializeLocalizer(localizer: SafeLocalizationDisambiguator): Unit = {
    writeState(TaskExecutor.LocalizerFile, localizer.toJson)
  }

  private def deserializeLocalizer(): SafeLocalizationDisambiguator = {
    SafeLocalizationDisambiguator.fromJson(readState(TaskExecutor.LocalizerFile))
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
    * Returns all schema types referenced in the inputs.
    */
  protected def getSchemas: Map[String, TSchema]

  /**
    * Evaluates input values and creates manifests for dxfuse and/or dxda to localize all
    * input files. Input files and directories are represented as URIs
    * (dx://proj-xxxx:file-yyyy::/A/B/C.txt). Input directories (VFolder and VListing types)
    * are updated recursively so that the listings are determined for all source directories.
    * Then, all the files are extracted from the inputs recursively and added to either the
    * dxfuse or dxda manifest, depending on the value of `streamFiles` and whether streaming
    * has been enabled individually for any of the inputs. The evaluated inputs and a mapping
    * from source URIs to local paths are stored to transient files for use by the next phase.
    * TODO: handle InitialWorkDir
    */
  def localizeInputs(): Unit = {
    if (logger.isVerbose) {
      trace(s"dxCompiler version: ${getVersion}")
      trace(s"Task source code:\n${jobMeta.sourceCode}", traceLengthLimit)
      printDirTree()
      trace(s"localizeInputs debugLevel=${logger.traceLevel}")
    }

    object PathsToLocalize {
      var virtualFiles: Vector[StringFileNode] = Vector.empty
      var localFiles: Set[LocalFileSource] = Set.empty
      var filesToStream: Set[AddressableFileNode] = Set.empty
      var filesToDownload: Set[AddressableFileNode] = Set.empty

      case class PathValueHandler(stream: Boolean) {
        // Updates all `VFolder` values (including those nested within `value`) to
        // set their `listing`s if missing. Also extracts all file `VFile` values into
        // a `PathsToLocalize` object as follows:
        // * If `contents` is non-empty, the file is added to `virtualFiles`
        // * If `address` resolves to a `LocalFileSource`, the file is added to `localFiles`
        // * Otherwise, if `stream` is true, the file is added to `filesToStream`, otherwise
        //   it is added to `filesToDownload`
        def updateListingsAndExtractFiles(value: Value,
                                          t: Option[Type] = None,
                                          optional: Boolean = false): Option[Value] = {
          // creates a listing of VFile/VFolders from a nested listing of FileSources
          def createListing(value: Vector[FileSource]): Vector[PathValue] = {
            value.map {
              case fs: AddressableFileSource if fs.isDirectory && fs.isListable =>
                val items = createListing(fs.listing())
                VFolder(fs.address, listing = Some(items))
              case fs: AddressableFileSource if fs.isDirectory =>
                VFolder(fs.address)
              case fs: LocalFileSource =>
                localFiles += fs
                VFile(fs.address)
              case fs: AddressableFileSource =>
                if (stream) {
                  filesToStream += fs
                } else {
                  filesToDownload += fs
                }
                VFile(fs.address)
            }
          }

          // updates a listing by ensuring all VFolders have their listing set
          def updateListing(items: Vector[PathValue]): Vector[PathValue] = {
            items.map { path =>
              updateListingsAndExtractFiles(path) match {
                case Some(p: PathValue) => p
                case other              => throw new Exception(s"expected PathValue, not ${other}")
              }
            }
          }

          def handleFileUri(uri: String): Unit = {
            fileResolver.resolve(uri) match {
              case fs if !fs.exists =>
                throw new Exception(s"File-type input does not exist: ${fs}")
              case fs if fs.isDirectory =>
                throw new Exception(s"File-type input is a directory: ${fs}")
              case fs: LocalFileSource => localFiles += fs
              case fs if stream        => filesToStream += fs
              case fs                  => filesToDownload += fs
            }
          }

          def getFolderListing(uri: String): Vector[FileSource] = {
            fileResolver.resolveDirectory(uri) match {
              case fs if !fs.exists =>
                throw new Exception(s"Directory-type input does not exist: ${fs}")
              case fs if !fs.isDirectory =>
                throw new Exception(s"Directory-type input is not a directory: ${fs}")
              case fs if !fs.isListable =>
                throw new Exception(s"Cannot get listing for Directory-type input: ${fs}")
              case fs => fs.listing(recursive = true)
            }
          }

          (t, value) match {
            case (Some(TFile) | None, file: VFile) if file.contents.nonEmpty =>
              virtualFiles :+= StringFileNode(file.contents.get, file.uri)
              Some(file.copy(contents = None))
            case (Some(TFile) | None, file: VFile) =>
              handleFileUri(file.uri)
              Some(file.copy(secondaryFiles = updateListing(file.secondaryFiles)))
            case (Some(TFile), VString(uri)) =>
              handleFileUri(uri)
              Some(VFile(uri))
            case (Some(TDirectory) | None, folder: VFolder) =>
              // fill in the listing if it is not provided
              Option
                .when(folder.listing.isEmpty)(getFolderListing(folder.uri))
                .filter(_.nonEmpty)
                .map { listing =>
                  val updatedListing = createListing(listing)
                  folder.copy(listing = Some(updatedListing))
                }
                .orElse(Some(folder))
            case (Some(TDirectory), VString(uri)) =>
              val listing = Some(getFolderListing(uri)).filter(_.nonEmpty).map(createListing)
              Some(VFolder(uri, listing = listing))
            case (Some(TDirectory) | None, listing: VListing) =>
              val updatedItems = updateListing(listing.items)
              Some(listing.copy(items = updatedItems))
            case _ => None
          }
        }
      }
    }

    // Make sure all VFolder inputs have their `listing`s set. We do this to ensure that
    // only the files present in the folder at the beginning of task execution are localized
    // into the inputs folder, and any changes to the directory during runtime aren't
    // synchronized to the worker. This step also extracts all the files from the inputs
    // that need to be localized.
    // TODO: it would be nice to extract dx:// links from VString values - this will happen
    //  in the case where the container is a dx file and being passed in as an input
    //  parameter - so that they could be downloaded using dxda. However, this would also
    //  require some way for the downloaded image tarball to be discovered and loaded. For now,
    //  we rely on DockerUtils to download the image (via DxFileSource, which uses the dx API
    //  to download the file).
    val inputs = getInputsWithDefaults.map {
      case (name, (t, v)) =>
        // whether to stream all the files/directories associated with this input
        val stream = streamFiles == StreamFiles.All ||
          (streamFiles == StreamFiles.PerFile && streamFileForInput(name))
        val handler = PathsToLocalize.PathValueHandler(stream)
        name -> (t, Value.transform(v, Some(t), handler.updateListingsAndExtractFiles))
    }
    serializeInputs(TaskExecutor.InitialInputsFile, getSchemas, inputs)

    // localize virtual files to /home/dnanexus/virtual
    logger.traceLimited(s"Virtual files = ${PathsToLocalize.virtualFiles}")
    val virtualUriToPath = PathsToLocalize.virtualFiles.map { fs =>
      val localPath = fs.localizeToDir(jobMeta.workerPaths.getVirtualFilesDir(ensureExists = true))
      fs.name -> localPath
    }.toMap

    // local files already have a path
    // TODO: should we use originalPath rather than address?
    val localUriToPath = PathsToLocalize.localFiles.map(fs => fs.address -> fs.canonicalPath).toMap

    // This object handles mapping FileSources to local paths and deals with file name
    // collisions in the manner specified by the WDL spec. All remote input files except those
    // streamed by dxfuse are placed in subfolders of the /home/dnanexus/inputs directory.
    val downloadLocalizer = SafeLocalizationDisambiguator.create(
        rootDir = jobMeta.workerPaths.getInputFilesDir(),
        separateDirsBySource = true,
        createDirs = true,
        disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
        logger = logger
    )
    // Build dxda manifest to localize all non-streaming remote files
    logger.traceLimited(s"Files to download: ${PathsToLocalize.filesToDownload}")
    val downloadFileSourceToPath: Map[AddressableFileNode, Path] =
      downloadLocalizer.getLocalPaths(PathsToLocalize.filesToDownload)
    val downloadUriToPath = downloadFileSourceToPath.map {
      case (fs, path) => fs.address -> path
    }
    DxdaManifestBuilder(dxApi, logger)
      .apply(downloadFileSourceToPath.collect {
        case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
      })
      .foreach {
        case DxdaManifest(manifestJs) =>
          FileUtils.writeFileContent(jobMeta.workerPaths.getDxdaManifestFile(),
                                     manifestJs.prettyPrint)
      }

    // build dxfuse manifest to localize all straming remote files and folders
    logger.traceLimited(s"Files to stream: ${PathsToLocalize.filesToStream}")
    // use a different localizer for the dxfuse mount point
    val streamingLocalizer = SafeLocalizationDisambiguator.create(
        rootDir = jobMeta.workerPaths.getDxfuseMountDir(),
        existingPaths = downloadLocalizer.getLocalizedPaths,
        separateDirsBySource = true,
        createDirs = false,
        disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
        logger = logger
    )
    val streamFileSourceToPath: Map[AddressableFileNode, Path] =
      streamingLocalizer.getLocalPaths(PathsToLocalize.filesToStream)
    val streamUriToPath = streamFileSourceToPath.map {
      case (fs, path) => fs.address -> path
    }
    DxfuseManifestBuilder(jobMeta.workerPaths, dxApi, logger)
      .apply(streamFileSourceToPath.collect {
        case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
      })
      .foreach {
        case DxfuseManifest(manifestJs) =>
          FileUtils.writeFileContent(jobMeta.workerPaths.getDxfuseManifestFile(),
                                     manifestJs.prettyPrint)
      }

    val allUriToPath = virtualUriToPath ++ localUriToPath ++ downloadUriToPath ++ streamUriToPath
    serializeUriToPath(allUriToPath)
    serializeLocalizer(downloadLocalizer)
  }

  /**
    * Finalizes all input files, folders, and listings.
    * For a file that is not in the context of a folder:
    * - If it is a local file, we link it into a disambiguation dir, naming it with
    *   its basename if it has one.
    * - Otherwise, if it has a basename, we create a link with the new name in
    *   the same folder to the original file, throwing an exception if there is
    *   a naming collision.
    * - If it has secondary files, they are linked into the same directory as
    *   the main file, creating any subfolders, and throwing an exception if there
    *   is a naming collision.
    * For a directory or listing, we create a new disambiguation dir and recursively
    * link in all the files it contains, creating any subfolders.
    */
  def finalizeInputs(): Unit = {
    val (schemas, inputs) = deserializeInputs(TaskExecutor.InitialInputsFile)
    val uriToPath = deserializeUriToPath()
    val localizer = deserializeLocalizer()
    logger.traceLimited(s"InstantiateCommand, env = ${inputs}, uriToPath = ${uriToPath}")

    object PathToUri {
      var pathToUri: Map[Path, String] = Map.empty
      val inputDirs =
        Vector(jobMeta.workerPaths.getInputFilesDir(), jobMeta.workerPaths.getDxfuseMountDir())

      def finalizeValue(value: Value, t: Option[Type], parent: Option[Path]): Option[Value] = {
        def finalizeListing(listing: Vector[PathValue], listingParent: Path): Vector[PathValue] = {
          listing.map { value =>
            finalizeValue(value, None, Some(listingParent)) match {
              case Some(p: PathValue) => p
              case other              => throw new Exception(s"expected PathValue, not ${other}")
            }
          }
        }

        (t, value) match {
          case (Some(TFile) | None, f: VFile) =>
            val path = uriToPath(f.uri)
            val parentDir = parent.getOrElse(path.getParent)
            val newPath =
              if (parent.nonEmpty || !inputDirs.exists(path.startsWith) || f.basename.nonEmpty) {
                // Either we are localizing into a specific directory, the path is a local or a
                // virtual file, or it is a remote file that needs its name changed - use the
                // localizer to determine the new path and then create a symbolic link.
                val sourcePath = f.basename.map(parentDir.resolve).getOrElse(path)
                val newPath = localizer.getLocalPath(fileResolver.fromPath(sourcePath))
                Files.createSymbolicLink(newPath, path)
              } else {
                path
              }
            pathToUri += (newPath -> f.uri)
            // recursively resolve the secondary files, linking them to be adjacent to the main file
            val newSecondaryFiles = finalizeListing(f.secondaryFiles, parentDir)
            Some(f.copy(uri = newPath.toString, secondaryFiles = newSecondaryFiles))
          case (Some(TDirectory) | None, f: VFolder) if f.listing.isDefined =>
            val fs = (fileResolver.resolveDirectory(f.uri), f.basename) match {
              case (fs, Some(basename)) =>
                fs.getParent
                  .map(_.resolve(basename))
                  .getOrElse(throw new Exception("cannot rename root directory"))
              case (fs, None) => fs
            }
            val path = localizer.getLocalPath(fs, parent)
            pathToUri += (path -> f.uri)
            val newListing = finalizeListing(f.listing.get, path)
            Some(f.copy(uri = path.toString, listing = Some(newListing)))
          case (Some(TDirectory) | None, l: VListing) =>
            // a listing is a virtual directory and is unique (i.e. we never want two
            // listings pointing at the same physical directory), so we source it from
            // a temp dir
            val tempDir = Files.createTempDirectory("listing")
            tempDir.toFile.deleteOnExit()
            val fs = fileResolver.fromPath(tempDir.resolve(l.basename), isDirectory = Some(true))
            val path = localizer.getLocalPath(fs, parent)
            val newListing = finalizeListing(l.items, path)
            Some(l.copy(items = newListing))
          case _ => None
        }
      }

      def handler(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
        finalizeValue(value, t, None)
      }
    }

    val finalizedInputs = inputs.map {
      case (name, (t, v)) =>
        name -> (t, Value.transform(v, Some(t), PathToUri.handler))
    }

    serializeInputs(TaskExecutor.LocalizedInputsFile, schemas, finalizedInputs)
    serializePathToUri(PathToUri.pathToUri)
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

  /**
    * Evaluates the command script and writes it to disk. Inputs are supplemented with
    * any local file paths created when evaluating the command script and are serialized
    * for use in the next phase.
    */
  def instantiateCommand(): Unit = {
    val (schemas, localizedInputs) = deserializeInputs(TaskExecutor.LocalizedInputsFile)
    logger.traceLimited(s"InstantiateCommand, env = ${localizedInputs}")
    // evaluate the command block and write the command script
    val updatedInputs = writeCommandScript(localizedInputs)
    // write the updated env to disk
    serializeInputs(TaskExecutor.LocalizedInputsWithCommandBindingsFile, schemas, updatedInputs)
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

  /**
    * Upload output files/directories and replace local paths with remote URIs
    * in the output values. An output file/folder may have been created on the
    * worker or it may be an input that originated either on the worker or from
    * the platform. We don't want to reupload files/directories that already
    * have an identical copy at the specified output location on the platform,
    * but otherwise we need to upload them.
    * TODO: handle output parameter tags and properties
    */
  def delocalizeOutputs(): Unit = {
    if (logger.isVerbose) {
      trace(s"delocalizeOutputs debugLevel=${logger.traceLevel}")
      printDirTree()
    }

    // load the cached localized inputs
    val (_, localizedInputs) = deserializeInputs(
        TaskExecutor.LocalizedInputsWithCommandBindingsFile
    )

    // evaluate output expressions
    val (localizedOutputs, tagsAndProperties) = evaluateOutputs(localizedInputs)

    // load the cached path mappings
    lazy val localPathToUri = deserializePathToUri()
    lazy val uriToSourcePath = deserializeUriToPath()
    lazy val dxInputDirs =
      Vector(jobMeta.workerPaths.getInputFilesDir(), jobMeta.workerPaths.getDxfuseMountDir())

    // does a local path correspond to a platform file?
    def isDxPath(path: Path): Boolean = {
      localPathToUri.get(path).exists { uri =>
        val sourcePath = uriToSourcePath(uri)
        dxInputDirs.exists(sourcePath.startsWith)
      }
    }

    // lazily build mapping of input path to VFolder
    lazy val localPathToInputFolder: Map[Path, VFolder] = {
      def handlePath(value: PathValue, folders: Map[Path, VFolder]): Map[Path, VFolder] = {
        value match {
          case f: VFolder =>
            fileResolver.resolveDirectory(f.uri) match {
              case local: LocalFileSource =>
                val foldersWithListing = f.listing
                  .map { listing =>
                    listing.foldLeft(folders) {
                      case (accu, path) => handlePath(path, accu)
                    }
                  }
                  .getOrElse(folders)
                foldersWithListing + (local.canonicalPath -> f)
              case _ => folders
            }
          case _ => folders
        }
      }

      def handleFolder(value: Value,
                       t: Option[Type],
                       optional: Boolean,
                       folders: Map[Path, VFolder]): Option[Map[Path, VFolder]] = {
        (t, value) match {
          case (Some(TDirectory) | None, f: VFolder) => Some(handlePath(f, folders))
          case _                                     => None
        }
      }

      localizedInputs.foldLeft(Map.empty[Path, VFolder]) {
        case (accu, (_, (t, v))) => Value.walk(v, Some(t), accu, handleFolder)
      }
    }

    // gets the listing from a folder input
    def getInputListing(localPath: Path): Option[Vector[PathValue]] = {
      localPathToInputFolder.get(localPath).flatMap(_.listing)
    }

    // compares the input and output listing for the same folder to see if they are the same
    def listingsEqual(localPath: Path, outputListing: Vector[PathValue]): Boolean = {
      // get all the paths in a single directory of a listing
      def getListingPaths(listing: Iterable[PathValue]): Set[Path] = {
        listing.toSet.flatMap {
          case f: VFile =>
            fileResolver.resolve(f.uri) match {
              case local: LocalFileSource =>
                Set(local.canonicalPath) ++ getListingPaths(f.secondaryFiles)
              case _ => Set.empty
            }
          case f: VFolder =>
            fileResolver.resolveDirectory(f.uri) match {
              case local: LocalFileSource => Set(local.canonicalPath)
              case _                      => Set.empty
            }
          case _ => Set.empty
        }
      }

      getInputListing(localPath).exists {
        case inputListing if inputListing.size != outputListing.size => false
        case inputListing =>
          val inputPaths = getListingPaths(inputListing)

          // An input and an output listing file/directory are equal if the resolve to the
          // same local path and the output item is not renamed (i.e. its basename is empty).
          // Directories are further required to have all their listing items be equal.
          def itemEqual(outputItem: PathValue): Boolean = {
            outputItem match {
              case f: VFile if f.basename.isDefined || f.contents.isDefined => false
              case f: VFile =>
                fileResolver.resolve(f.uri) match {
                  case local: LocalFileSource =>
                    inputPaths.contains(local.canonicalPath) && f.secondaryFiles.forall(itemEqual)
                  case _ => false
                }
              case f: VFolder if f.basename.isDefined => false
              case f: VFolder =>
                fileResolver.resolveDirectory(f.uri) match {
                  case local: LocalFileSource =>
                    inputPaths.contains(local.canonicalPath) && f.listing.forall(listing =>
                      listingsEqual(local.canonicalPath, listing)
                    )
                  case _ => false
                }
            }
          }

          outputListing.forall(itemEqual)
      }
    }

    // creates a listing of VFile/VFolders from a nested listing of FileSources
    def createListing(value: Vector[FileSource]): Vector[PathValue] = {
      value.map {
        case fs: LocalFileSource if fs.isDirectory =>
          VFolder(fs.address, listing = Some(createListing(fs.listing())))
        case fs: LocalFileSource => VFile(fs.address)
        case other               => throw new Exception(s"unexpected output listing item ${other}")
      }
    }

    // update outputs - make sure all folders have listings
    def setFolderListings(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      (t, value) match {
        case (Some(TDirectory) | None, f: VFolder) if f.listing.isEmpty =>
          fileResolver.resolveDirectory(f.uri) match {
            case local: LocalFileSource =>
              Some(f.copy(listing = getInputListing(local.canonicalPath).orElse {
                if (!Files.exists(local.canonicalPath)) {
                  throw new Exception(
                      s"Directory-type output does not exist: ${local.canonicalPath}"
                  )
                }
                Some(createListing(local.listing(recursive = true)))
              }))
            case _ => None
          }
        case _ => None
      }
    }

    val updatedOutputs = localizedOutputs.map {
      case (name, (t, v)) => name -> (t, Value.transform(v, Some(t), setFolderListings))
    }

    val defaultDestParent = if (jobMeta.useManifests) {
      s"${jobMeta.manifestFolder}/"
    } else {
      "/"
    }

    case class Uploads(tags: Set[String], properties: Map[String, String]) {
      private var files: Vector[FileUpload] = Vector.empty
      private var virtualFiles: Map[String, Path] = Map.empty
      private var directories: Map[Path, String] = Map.empty

      def addFile(localPath: Path, destination: String): Unit = {
        files :+= FileUpload(source = localPath,
                             destination = Some(destination),
                             tags = tags,
                             properties = properties)
      }

      def addVirtualFile(content: String, uri: String, destination: String): Unit = {
        val localPath = jobMeta.workerPaths.getVirtualFilesDir(ensureExists = true).resolve(uri)
        val canonicalPath = FileUtils.writeFileContent(localPath, content)
        virtualFiles += (uri -> canonicalPath)
        addFile(canonicalPath, destination)
      }

      def addDirectory(localPath: Path, destination: String): Unit = {
        directories += (localPath -> destination)
      }

      def asTuple: (Vector[FileUpload], Map[String, Path], Map[Path, String]) =
        (files, virtualFiles, directories)
    }

    // extract all files to upload
    def extractPaths(value: Value,
                     t: Option[Type],
                     optional: Boolean,
                     ctx: Uploads): Option[Uploads] = {
      def handlePath(innerValue: PathValue,
                     innerType: Option[Type],
                     parent: Option[String] = None): Unit = {
        (innerType, innerValue) match {
          case (Some(TFile) | None, f: VFile) if f.contents.isDefined =>
            // a file literal
            val name = f.basename.getOrElse(f.uri)
            val destFolder = parent.getOrElse(defaultDestParent)
            ctx.addVirtualFile(f.contents.get, f.uri, s"${destFolder}${name}")
            f.secondaryFiles.foreach { path =>
              handlePath(path, None, Some(destFolder))
            }
          case (Some(TFile) | None, f: VFile) =>
            fileResolver.resolve(f.uri) match {
              case local: LocalFileSource
                  if parent.isDefined || f.basename.isDefined || !isDxPath(local.canonicalPath) =>
                // Either the file was created by the command script, it originated as
                // a local or virtual file, it needs to be uploaded with a new name, or
                // we are uploading a listing/directory and so need to make a new copy
                // of the file regardless of its origin
                val name = f.basename.getOrElse(local.canonicalPath.getFileName.toString)
                val destFolder = parent.getOrElse(defaultDestParent)
                ctx.addFile(local.canonicalPath, s"${destFolder}${name}")
                f.secondaryFiles.foreach { path =>
                  handlePath(path, None, Some(destFolder))
                }
              case _ => ()
            }
          case (Some(TDirectory) | None, f: VFolder) =>
            fileResolver.resolveDirectory(f.uri) match {
              case local: LocalFileSource
                  if Files.exists(local.canonicalPath) && (
                      parent.isDefined ||
                        f.basename.isDefined ||
                        !isDxPath(local.canonicalPath) ||
                        f.listing.exists(listing => !listingsEqual(local.canonicalPath, listing))
                  ) =>
                // Either the directory was created by the command script, it originated as
                // a local directory, it needs to be uploaded with a new name, its listing
                // is different, or we are uploading a subdirectory of a listing/directory
                // and so need to make a new copy of the this directory regardless of its origin
                val name = f.basename.getOrElse(local.canonicalPath.getFileName.toString)
                val destFolder = s"${parent.getOrElse(defaultDestParent)}${name}/"
                ctx.addDirectory(local.canonicalPath, destFolder)
                f.listing.get.foreach { item =>
                  handlePath(item, None, Some(destFolder))
                }
              case _ => ()
            }
          case (Some(TDirectory) | None, l: VListing) =>
            val destFolder = s"${parent.getOrElse(defaultDestParent)}${l.basename}/"
            l.items.foreach { item =>
              handlePath(item, None, Some(destFolder))
            }
          case _ =>
            throw new Exception(s"Invalid path value ${innerValue} for type ${innerType}")
        }
      }
      (t, value) match {
        case (Some(TFile) | Some(TDirectory) | None, p: PathValue) =>
          handlePath(p, t)
          Some(ctx)
        case _ => None
      }
    }

    // upload files
    val (filesToUpload, virtualFiles, directories) = updatedOutputs.map {
      case (name, (t, v)) =>
        val (tags, properties) =
          tagsAndProperties.getOrElse(name, (Set.empty[String], Map.empty[String, String]))
        Value.walk(v, Some(t), Uploads(tags, properties), extractPaths).asTuple
    }.unzip3
    val uploadedFiles = filesToUpload.flatten match {
      case itr if itr.nonEmpty =>
        dxApi.uploadFiles(files = itr,
                          waitOnUpload = waitOnUpload,
                          maxConcurrent = TaskExecutor.MaxConcurrentUploads)
      case _ => Map.empty[Path, DxFile]
    }

    val delocalizedOutputs = if (uploadedFiles.nonEmpty) {
      val uploadedVirtualFiles = virtualFiles.flatten.toMap
      val uploadedDirectories = directories.flatten.toMap
      // For folders/listings, we have to determine the project output folder,
      // because the platform won't automatically update directory URIs like
      // it does for data object links. The output folder is either the manifest
      // folder or the concatination of the job/analysis output folder and the
      // container output folder.
      // TODO: currently, this only works for one level of job nesting -
      //  we probably need to get all the output folders all the way up
      //  the job tree and join them all
      val defaultFolderParent = jobMeta.projectOutputFolder match {
        case folder if folder.endsWith("/") => folder
        case folder                         => s"${folder}/"
      }

      // replace local paths with remote URIs
      def handlePath(value: PathValue,
                     t: Option[Type] = None,
                     optional: Boolean = false,
                     parent: Option[String] = None): PathValue = {
        case (Some(TFile) | None, f: VFile) if f.contents.isDefined =>
          // a virtual file that was uploaded - replace the URI with the dx URI and
          // unset the contents
          val dxFile = uploadedFiles(uploadedVirtualFiles(f.uri))
          VFile(dxFile.asUri, secondaryFiles = f.secondaryFiles.map(handlePath(_)))
        case (Some(TFile) | None, f: VFile) =>
          val newSecondaryFiles = f.secondaryFiles.map(handlePath(_))
          fileResolver.resolve(f.uri) match {
            case local: LocalFileSource if optional && !Files.exists(local.canonicalPath) =>
              VNull
            case local: LocalFileSource if uploadedFiles.contains(local.canonicalPath) =>
              // an output file that was uploaded
              f.copy(uri = uploadedFiles(local.canonicalPath).asUri,
                     secondaryFiles = newSecondaryFiles)
            case local: LocalFileSource if localPathToUri.contains(local.canonicalPath) =>
              // an output file that was an input file
              f.copy(uri = localPathToUri(local.canonicalPath), secondaryFiles = newSecondaryFiles)
            case _: LocalFileSource =>
              throw new Exception(s"file was not delocalized: ${f}")
            case _ => f
          }
        case (Some(TDirectory) | None, f: VFolder) =>
          fileResolver.resolveDirectory(f.uri) match {
            case local: LocalFileSource if optional && !Files.exists(local.canonicalPath) =>
              VNull
            case local: LocalFileSource if uploadedDirectories.contains(local.canonicalPath) =>
              // a folder that was uploaded
              val folderParent = parent.getOrElse(defaultFolderParent)
              val name = f.basename.getOrElse(local.canonicalPath.getFileName.toString)
              val folder = s"${folderParent}${name}"
              VFolder(folder, listing = f.listing.map(_.map(handlePath(_, parent = Some(folder)))))
            case local: LocalFileSource
                if localPathToUri.contains(local.canonicalPath) && parent.isEmpty =>
              // an output folder that was an input folder
              VFolder(localPathToUri(local.canonicalPath),
                      listing = getInputListing(local.canonicalPath))
            case _: LocalFileSource =>
              throw new Exception(s"folder was not delocalized: ${f}")
            case _ => f
          }
        case (Some(TDirectory) | None, l: VListing) =>
          val folderParent = parent.getOrElse(defaultFolderParent)
          val folder = s"${folderParent}${l.basename}"
          VFolder(folder, listing = Some(l.items.map(handlePath(_, parent = Some(folder)))))
        case _ =>
          throw new Exception(s"Could not delocalize ${value} as ${t}")
      }

      def delocalizePaths(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
        (t, value) match {
          case (Some(TFile) | Some(TDirectory) | None, p: PathValue) =>
            Some(handlePath(p, t, optional))
          case _ => None
        }
      }

      updatedOutputs.map {
        case (name, (t, v)) =>
          name -> (t, Value.transform(v, Some(t), delocalizePaths))
      }
    } else {
      updatedOutputs
    }

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
          case TaskAction.LocalizeInputs     => localizeInputs()
          case TaskAction.FinalizeInputs     => finalizeInputs()
          case TaskAction.InstantiateCommand => instantiateCommand()
          case TaskAction.DelocalizeOutputs  => delocalizeOutputs()
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
