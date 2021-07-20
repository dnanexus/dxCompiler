package dx.executor

import java.nio.file.{Files, Path, Paths}
import dx.api.{DxFile, DxJob, InstanceTypeRequest}
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
import dx.util.protocols.{DxFileAccessProtocol, DxFileSource, DxFolderSource}
import spray.json._
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

object TaskAction extends Enum {
  type TaskAction = Value
  val CheckInstanceType, LocalizeInputs, FinalizeInputs, InstantiateCommand, DelocalizeOutputs,
      Relaunch = Value
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

  private def serializeInputs(schemas: Map[String, Type],
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
    val json = JsObject(
        "schemas" -> JsObject(newSchemasJs),
        "localizedInputs" -> JsObject(inputsJs)
    )
    FileUtils.writeFileContent(
        jobMeta.workerPaths.getSerializedInputsFile(ensureParentExists = true),
        json.prettyPrint
    )
  }

  private def deserializeInputs(): (Map[String, Type], Map[String, (Type, Value)]) = {
    val localizedInputsFile = jobMeta.workerPaths.getSerializedInputsFile()
    val (schemasJs, inputsJs) = FileUtils
      .readFileContent(localizedInputsFile)
      .parseJson match {
      case env: JsObject =>
        env
          .getFields("schemas", "localizedInputs") match {
          case Seq(JsObject(schemas), JsObject(inputs)) => (schemas, inputs)
          case _ =>
            throw new Exception(s"Malformed environment serialized to disk ${localizedInputsFile}")
        }
      case _ =>
        throw new Exception(s"Malformed environment serialized to disk ${localizedInputsFile}")
    }
    val schemas = TypeSerde.deserializeSchemas(schemasJs)
    inputsJs
      .foldLeft((schemas, Map.empty[String, (Type, Value)])) {
        case ((schemaAccu, paramAccu), (key, obj: JsObject)) =>
          obj.getFields("type", "value") match {
            case Seq(typeJs, valueJs) =>
              val (irType, newSchemas) = TypeSerde.deserialize(typeJs, schemaAccu)
              val irValue = ValueSerde.deserializeWithType(valueJs, irType, key)
              (newSchemas, paramAccu + (key -> (irType, irValue)))
            case _ =>
              throw new Exception(
                  s"Malformed environment serialized to disk ${localizedInputsFile}"
              )
          }
        case _ =>
          throw new Exception(s"Malformed environment serialized to disk ${localizedInputsFile}")
      }
  }

  private def serializeUriToPath(uriToPath: Map[String, Path]): Unit = {
    val utiToPathJs: Map[String, JsValue] = uriToPath.map {
      case (uri, path) => uri -> JsString(path.toString)
    }
    FileUtils.writeFileContent(
        jobMeta.workerPaths.getSerializedUriToPathFile(ensureParentExists = true),
        JsObject(utiToPathJs).prettyPrint
    )
  }

  private def deserializeUriToPath(): Map[String, Path] = {
    val uriToPathFile = jobMeta.workerPaths.getSerializedUriToPathFile()
    FileUtils.readFileContent(uriToPathFile).parseJson match {
      case JsObject(uriToPath) =>
        uriToPath.map {
          case (uri, JsString(path)) => uri -> Paths.get(path)
          case other                 => throw new Exception(s"unexpected path ${other}")
        }
      case _ =>
        throw new Exception(s"Malformed environment serialized to disk ${uriToPathFile}")
    }
  }

  private def serializePathToUri(pathToUri: Map[Path, String]): Unit = {
    val pathToUriJs: Map[String, JsValue] = pathToUri.map {
      case (path, uri) => path.toString -> JsString(uri)
    }
    FileUtils.writeFileContent(
        jobMeta.workerPaths.getSerializedPathToUriFile(ensureParentExists = true),
        JsObject(pathToUriJs).prettyPrint
    )
  }

  private def deserializePathToUri(): Map[Path, String] = {
    val pathToUriFile = jobMeta.workerPaths.getSerializedPathToUriFile()
    FileUtils.readFileContent(pathToUriFile).parseJson match {
      case JsObject(pathToUri) =>
        pathToUri.map {
          case (path, JsString(uri)) => Paths.get(path) -> uri
          case other                 => throw new Exception(s"unexpected uri ${other}")
        }
      case _ =>
        throw new Exception(s"Malformed environment serialized to disk ${pathToUriFile}")
    }
  }

  private def serializeLocalizer(localizer: SafeLocalizationDisambiguator): Unit = {
    FileUtils.writeFileContent(
        jobMeta.workerPaths.getSerializedLocalizerFile(ensureParentExists = true),
        localizer.toJson.prettyPrint
    )
  }

  private def deserializeLocalizer(): SafeLocalizationDisambiguator = {
    SafeLocalizationDisambiguator.fromJson(
        FileUtils.readFileContent(jobMeta.workerPaths.getSerializedLocalizerFile()).parseJson
    )
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

          (t, value) match {
            case (Some(TFile) | None, file: VFile) if file.contents.nonEmpty =>
              virtualFiles :+= StringFileNode(file.contents.get, file.uri)
              Some(file.copy(contents = None))
            case (Some(TFile) | None, file: VFile) =>
              fileResolver.resolve(file.uri) match {
                case fs if !fs.exists =>
                  throw new Exception(s"File-type input does not exist: ${fs}")
                case fs if fs.isDirectory =>
                  throw new Exception(s"File-type input is a directory: ${fs}")
                case fs: LocalFileSource => localFiles += fs
                case fs if stream        => filesToStream += fs
                case fs                  => filesToDownload += fs
              }
              val updatedSecondaryFiles = updateListing(file.secondaryFiles)
              Some(file.copy(secondaryFiles = updatedSecondaryFiles))
            case (Some(TDirectory) | None, folder: VFolder) =>
              // fill in the listing if it is not provided
              Option
                .when(folder.listing.isEmpty)(fileResolver.resolveDirectory(folder.uri))
                .map {
                  case fs if !fs.exists =>
                    throw new Exception(s"Directory-type input does not exist: ${fs}")
                  case fs if !fs.isDirectory =>
                    throw new Exception(s"Directory-type input is not a directory: ${fs}")
                  case fs if !fs.isListable =>
                    throw new Exception(s"Cannot get listing for Directory-type input: ${fs}")
                  case fs => fs.listing(recursive = true)
                }
                .filter(_.nonEmpty)
                .map { listing =>
                  // fill in the listing if it is not provided
                  val updatedListing = createListing(listing)
                  folder.copy(listing = Some(updatedListing))
                }
                .orElse(Some(folder))
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
    serializeInputs(getSchemas, inputs)

    // localize virtual files to /home/dnanexus/virtual
    logger.traceLimited(s"Virtual files = ${PathsToLocalize.virtualFiles}")
    val virtualUriToPath = PathsToLocalize.virtualFiles.map { fs =>
      val localPath = fs.localizeToDir(jobMeta.workerPaths.getVirtualFilesDir(ensureExists = true))
      fs.name -> localPath
    }.toMap

    // local files already have a path
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

  //  val uriToPath: Map[String, String] = (fileSourceToPath ++ folderSourceToPath).map {
  //    case (dxFs: DxFileSource, path)       => dxFs.address -> path.toString
  //    case (dxFs: DxFolderSource, path)     => dxFs.address -> path.toString
  //    case (localFs: LocalFileSource, path) => localFs.originalPath.toString -> path.toString
  //    case (other, _)                       => throw new RuntimeException(s"unsupported file source ${other}")
  //  }.toMap

  def finalizeInputs(): Unit = {
    val (schemas, inputs) = deserializeInputs()
    val uriToPath = deserializeUriToPath()
    val localizer = deserializeLocalizer()
    logger.traceLimited(s"InstantiateCommand, env = ${inputs}, uriToPath = ${uriToPath}")

    // Go through all the inputs and finalize files, folders and listings.
    // For a file that is not in the context of a folder:
    // * If it is a local file, we link it into a disambiguation dir, naming it with
    //   its basename if it has one.
    // * Otherwise, if it has a basename, we create a link with the new name in
    //   the same folder to the original file, throwing an exception if there is
    //   a naming collision.
    // * If it has secondary files, they are linked into the same directory as
    //   the main file, creating any subfolders, and throwing an exception if there
    //   is a naming collision.
    // For a directory or listing, we create a new disambiguation dir and recursively
    // link in all the files it contains, creating any subfolders.
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

    serializeInputs(schemas, finalizedInputs)
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
    * Updates inputs with localized paths, then evaluates the command script and
    * writes it to disk. Inputs are supplemented with any local file paths created
    * when evaluating the command script and are serialized for use in the next phase.
    */
  def instantiateCommand(): Unit = {
    val (schemas, localizedInputs) = deserializeInputs()
    logger.traceLimited(s"InstantiateCommand, env = ${localizedInputs}")
    // evaluate the command block and write the command script
    val updatedInputs = writeCommandScript(localizedInputs)
    // write the updated env to disk
    serializeInputs(schemas, updatedInputs)
  }

  /**
    * Evaluates the outputs of the task. Returns mapping of output parameter
    * name to (type, value, Set of tags, and Map of properties), where tags
    * and properites only apply to output files (or collections thereof).
    * @param localizedInputs the job inputs, with files localized to the worker
    */
  protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): (Map[String, (Type, Value)], Map[String, (Set[String], Map[String, String])])

  def delocalizeOutputs(): Unit = {
    if (logger.isVerbose) {
      trace(s"delocalizeOutputs debugLevel=${logger.traceLevel}")
      printDirTree()
    }

    // load the cached localized inputs and path mappings
    val (_, localizedInputs) = deserializeInputs()
    val localPathToUri = deserializePathToUri()
    val uriToSourcePath = deserializeUriToPath()

    // evaluate output expressions
    val (localizedOutputs, tagsAndProperties) = evaluateOutputs(localizedInputs)

    // extract all files to upload
    def extractPaths(value: Value,
                     t: Option[Type],
                     optional: Boolean,
                     ctx: Vector[Path]): Option[Vector[Path]] = {
      (t, value) match {
        case (Some(TFile) | None, f: VFile) if f.contents.isDefined =>
          // a file literal
          val name = f.basename.getOrElse(f.uri)

        case _ => None
      }
    }

    val outputFiles = localizedOutputs.map {
      case (name, (t, v)) => name -> Value.walk(v, Some(t), Vector.empty[Path], extractPaths)
    }

    // upload files
    jobMeta.uploadOutputFiles(outputFiles, tagsAndProperties)

    // replace local paths with remote URIs

  }

  /**
    * Collection of files associated with an input value.
    * @param virutalFiles VFiles that have a `contents` attribute - these are
    *                     localized by writing the contents to a file
    * @param realFiles local/remote files
    * @param folders local/remove folders - these are localized by streaming/
    *                downloading the folder recursively
    */
  private case class PathValues(virutalFiles: Vector[VFile] = Vector.empty,
                                realFiles: Vector[(AddressableFileNode, VFile)] = Vector.empty,
                                folders: Vector[(AddressableFileSource, VFolder)] = Vector.empty)

  private lazy val virtualOutputDir =
    Files.createTempDirectory(jobMeta.workerPaths.getOutputFilesDir(ensureExists = true), "virtual")

  private def extractPaths(v: Value, fillInListings: Boolean): PathValues = {
    def extractArray(values: Vector[Value]): PathValues = {
      values.map(extractPaths(_, fillInListings)).foldLeft(PathValues()) {
        case (pathsAccu, paths) =>
          PathValues(pathsAccu.virutalFiles ++ paths.virutalFiles,
                     pathsAccu.realFiles ++ paths.realFiles,
                     pathsAccu.folders ++ paths.folders)
      }
    }
    def createListing(value: Vector[FileSource]): Vector[PathValue] = {
      value.collect {
        case fs: AddressableFileSource if fs.isDirectory && fs.isListable =>
          VFolder(fs.address, listing = Some(createListing(fs.listing())))
        case fs: AddressableFileSource if fs.isDirectory => VFolder(fs.address)
        case fs: AddressableFileSource                   => VFile(fs.address)
      }
    }
    v match {
      case file: VFile if file.contents.nonEmpty =>
        PathValues(virutalFiles = Vector(file))
      case file: VFile =>
        val secondaryPaths = extractArray(file.secondaryFiles)
        val fileSource = fileResolver.resolve(file.uri)
        secondaryPaths.copy(realFiles = secondaryPaths.realFiles :+ (fileSource -> file))
      case folder: VFolder =>
        val folderSource = fileResolver.resolveDirectory(folder.uri)
        // fill in the listing if it is not provided
        val updatedFolder = if (folder.listing.isEmpty && folderSource.isListable) {
          folder.copy(listing = Some(createListing(folderSource.listing(recursive = true))))
        } else {
          folder
        }
        PathValues(folders = Vector(folderSource -> updatedFolder))
      case listing: VListing => extractArray(listing.items)
      case VArray(items)     => extractArray(items)
      case VHash(m)          => extractArray(m.values.toVector)
      case _                 => PathValues()
    }
  }

  // TODO: handle secondary files - put in same output folder as the primary
  private def extractOutputPaths(
      name: String,
      v: Value,
      t: Type
  ): (Vector[(AddressableFileNode, VFile)], Vector[(AddressableFileSource, VFolder)]) = {
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
    ): (Vector[(LocalFileSource, VFile)], Vector[(LocalFileSource, VFolder)]) = {
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
          (Vector(fileResolver.fromPath(resolved) -> f), Vector.empty)
        case (TFile, f: VFile) =>
          (getLocalFileSource(innerName, fileResolver.resolve(f.uri), optional)
             .map(_ -> f)
             .toVector,
           Vector.empty)
        case (TFile, VString(path)) =>
          (getLocalFileSource(innerName, fileResolver.resolve(path), optional)
             .map(_ -> VFile(path))
             .toVector,
           Vector.empty)
        case (TDirectory, f: VFolder) =>
          (Vector.empty,
           getLocalFileSource(innerName, fileResolver.resolveDirectory(f.uri), optional)
             .map(_ -> f)
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
                case (fs, file) =>
                  getLocalFileSource(s"${innerName}.${key}", fs, optional = true)
                    .map(_ -> file)
                    .toVector
              }
              val folders = fileSources.folders.flatMap {
                case (fs, folder) =>
                  getLocalFileSource(s"${innerName}.${key}", fs, optional = true)
                    .map(_ -> folder)
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

  def delocalizeOutputs(): Unit = {
    if (logger.isVerbose) {
      trace(s"Epilog debugLevel=${logger.traceLevel}")
      printDirTree()
    }

    // load the cached localized inputs - folders have no listings
    val (_, localizedInputs) = deserializeInputs()
    // evaluate output expressions
    val localizedOutputs = evaluateOutputs(localizedInputs)
    // extract files/folders from the outputs
    val (localOutputFiles, localOutputFolders) = localizedOutputs.map {
      case (name, (irType, irValue)) => extractOutputPaths(name, irValue, irType)
    }.unzip

    <<<<<<< HEAD
    // Build a map of all the addresses in the output values that might
    // map to the same (absolute) local path. An output may be a file/folder
    // that was an input (or nested within an input folder) - these do not
    // need to be re-uploaded. Since the input directories are read-only (i.e.
    // there is no possibility the user created a new file there), we just need
    // to check whether the local file/folder has a path starting with one of
    // the input directories to determine if it needs to be uploaded.
    val (fileSourceToPath, folderSourceToPath) = deserializeFileSourceToPath()
    val inputDirs =
      Set(jobMeta.workerPaths.getInputFilesDir(), jobMeta.workerPaths.getDxfuseMountDir())
    val delocalizingPathToFile = localOutputFiles.flatten.toMap.collect {
      case (local: LocalFileSource, file) if !inputDirs.exists(local.canonicalPath.startsWith) =>
        local.canonicalPath -> file
    }
    val delocalizingPathToFolder = localOutputFolders.flatten.toMap.collect {
      case (local: LocalFileSource, folder) if !inputDirs.exists(local.canonicalPath.startsWith) =>
        local.canonicalPath -> folder
    }

    // Upload the files/directories, and map their local paths to their remote URIs.
    // If using manifests, we need to upload the files directly to the project.
    val destFolder = Option.when(jobMeta.useManifests)(jobMeta.manifestFolder)

    def ensureAbsolute(basename: String): String = {
      destFolder.map(f => s"${f}/${basename}").getOrElse(s"/${basename}")
    }

    // Gets all the file/folder paths from a listing
    // TODO: deal with secondary files
    // TODO: deal with basenames
    def getListingPaths(root: Path, pathValue: PathValue): Set[Path] = {
      def getPath(uri: String): Path = {
        // uri should be a relative path or an absolute path that is a
        // child of root
        fileResolver.resolve(uri) match {
          case local: LocalFileSource if local.canonicalPath.isAbsolute =>
            local.canonicalPath
          case local: LocalFileSource =>
            root.resolve(local.canonicalPath)
          case other =>
            throw new Exception(s"not a LocalFileSource ${other}")
        }
      }
      pathValue match {
        case f: VFile if f.contents.isEmpty => Set(getPath(f.uri))
        case f: VFolder if f.listing.isDefined =>
          f.listing.get.flatMap(l => getListingPaths(getPath(f.uri), l)).toSet
        case f: VFolder => Set(getPath(f.uri))
        case other =>
          throw new Exception(s"unsupported listing item ${other}")
      }
    }

    case class Directory(path: Path,
                         var files: Map[String, DxFile] = Map.empty,
                         var subdirs: Map[String, Directory] = Map.empty) {
      def addFile(name: String, dxFile: DxFile): Unit = {
        files += (name -> dxFile)
      }

      def addSubdir(name: String, path: Path): Directory = {
        val subdir = Directory(path)
        subdirs += (name -> subdir)
        subdir
      }
    }

    // creates a hierarchical listing from a collection of DxFiles with no
    // guaranteed ordering
    def createListing(root: Path, files: Map[Path, DxFile]): Directory = {
      val rootDir = Directory(root)
      var dirs: Map[Path, Directory] = Map(root -> rootDir)

      def getDir(parent: Path): Directory = {
        if (parent == null) {
          throw new Exception(s"${parent} is not a child of ${root}")
        }
        dirs.getOrElse(parent, {
          val grandparent = getDir(parent)
          val dir = grandparent.addSubdir(parent.getFileName.toString, parent)
          dirs += (parent -> dir)
          dir
        })
      }

      files.foreach {
        case (path, dxFile) => getDir(path.getParent).addFile(path.getFileName.toString, dxFile)
      }

      rootDir
    }

    val dxProtocol = fileResolver.getProtocolForScheme("dx") match {
      case proto: DxFileAccessProtocol => proto
      case other =>
        throw new Exception(s"unexpected dx protocol ${other}")
    }

    // upload files
    val fileToDelocalizedUri: Map[Path, DxFileSource] = fileUploader
      .uploadFilesWithDestination(
          delocalizingPathToFile.map {
            case (path, file) =>
              val basename = file.basename.getOrElse(path.getFileName.toString)
              path -> ensureAbsolute(basename)
          },
          wait = waitOnUpload
      )
      .map {
        case (path, dxFile) => path -> dxProtocol.fromDxFile(dxFile)
      }

    // upload folders - use listing to restrict the files/subfolders that are uploaded
    val folderToDelocalizedUriAndListing: Map[Path, (DxFolderSource, Directory)] = {
      val (dirs, listings) = delocalizingPathToFolder.map {
        case (path, folder) =>
          val basename = folder.basename.getOrElse(path.getFileName.toString)
          // if the VFolder has a listing, restrict the upload to only those files/folders
          val listingPaths = getListingPaths(path, folder)
          (path -> ensureAbsolute(basename),
           Option.unless(listingPaths.isEmpty)(path -> listingPaths))
      }.unzip
      val jobFolder = Paths.get(jobMeta.folder)
      fileUploader
        .uploadDirectoriesWithDestination(dirs.toMap,
                                          wait = waitOnUpload,
                                          listings = listings.flatten.toMap)
        .map {
          case (path, (projectId, folder, files)) =>
            // the folder is in the job/analysis container, which is relative to the
            // job/analysis output folder, so we need to resolve here
            // TODO: currently, this only works for one level of job nesting -
            //  we probably need to get all the output folders all the way up
            //  the job tree and join them all
            val containerFolder = Paths.get(folder)
            val baseFolder = jobFolder.resolve(if (containerFolder.isAbsolute) {
              Paths.get("/").relativize(containerFolder)
            } else {
              containerFolder
            })
            // create the output listing - this is necessary because the output
            // object must contain dx links for all the uploaded files or the
            // files will not be added to the job output (and thus not cloned
            // into the project)
            path -> (
                dxProtocol.fromDxFolder(projectId, baseFolder.toString),
                createListing(path, files)
            )
        }
    }
    =======
    // extract files from the outputs
    val localOutputFileSources = localizedOutputs.flatMap {
      case (name, (irType, irValue, tags, properties)) =>
        extractOutputFiles(name, irValue, irType).map(fs => (fs, tags, properties))
    }.toVector

    // Build a map of all the string values in the output values that might
    // map to the same (absolute) local path. Some of the outputs may be files
    // that were inputs (in `fileSourceToPath`) - these do not need to be
    // re-uploaded. The `localPath`s will be the same but the `originalPath`s
    // may be different.
    val inputAddresses = fileSourceToPath.keySet.map(_.address)
    val inputPaths = fileSourceToPath.values.toSet
    val filesToUpload = localOutputFileSources.collect {
      case (local: LocalFileSource, tags, properties)
          if !(
              inputAddresses.contains(local.address) || inputPaths.contains(local.canonicalPath)
          ) =>
        // if using manifests, we need to upload the files directly to the project
        val dest = if (jobMeta.useManifests) {
          Some(s"${jobMeta.manifestFolder}/${local.canonicalPath.getFileName.toString}")
        } else {
          None
        }
        local.address -> FileUpload(local.canonicalPath, dest, tags, properties)
    }.toMap

    val delocalizingValueToPath = filesToUpload.map {
      case (addr, FileUpload(path, _, _, _)) => addr -> path
    }

    // upload the files, and map their local paths to their remote URIs
    val delocalizedPathToUri =
      fileUploader.upload(filesToUpload.values.toSet, wait = waitOnUpload).map {
        case (path, dxFile) => path -> dxFile.asUri
      }
    >>>>>>> develop

    // Replace the local paths in the output values with URIs. For a file/folder
    // that is located in the inputs folder, we search in the localized inputs
    // for the file/folder or its nearest ancestor directory. For files/folders
    // that were generated on the worker, we need to do two look-ups: first to
    // get the absolute Path associated with the file value (which may be relative
    // or absolute), and second to get the URI associated with the Path. Returns
    // an Optional[String] because optional outputs may not be specified.
    val inputFileToSource = fileSourceToPath
      .collect {
        case (fs: AddressableFileNode, path) =>
          Map(fs.address -> fs, path.toString -> fs)
      }
      .flatten
      .toMap
    val delocalizingFileUriToPath = delocalizingPathToFile.map {
      case (path, file) => file.uri -> path
    }
    val inputFolderToSource = folderSourceToPath
      .collect {
        case (fs: AddressableFileSource, path) =>
          Map(fs.address -> fs, path.toString -> fs)
      }
      .flatten
      .toMap
    val delocalizingFolderUriToPath = delocalizingPathToFolder.map {
      case (path, folder) => folder.uri -> path
    }

    // Replace the URIs with remote file paths
    def pathTranslator(v: Value, t: Option[Type], optional: Boolean): Option[PathValue] = {
      def resolveInputSource(fs: AddressableFileSource): Option[AddressableFileSource] = {
        fs.getParent.flatMap { parent =>
          inputFolderToSource
            .get(parent.address)
            .orElse(resolveInputSource(parent))
            .map(_.resolve(fs.name))
        }
      }
      def translateFileUri(uri: String): Option[String] = {
        inputFileToSource
          .get(uri) // look for the local URI directly in the inputs
          .orElse(resolveInputSource(fileResolver.resolve(uri))) // look for the local URI nested in an input folder
          .orElse( // look for the local URI in the set that were uploaded
              delocalizingFileUriToPath
                .get(uri) match {
                case Some(path) => fileToDelocalizedUri.get(path)
                case _          => None
              }
          ) match {
          case Some(fs)         => Some(fs.address)
          case None if optional => None
          case _                => throw new Exception(s"Did not localize file ${uri}")
        }
      }
      def translateFolderUri(uri: String): Option[(String, Directory)] = {
        inputFolderToUri
          .get(uri)
          .orElse(delocalizingFolderUriToPath.get(uri) match {
            case Some(path) => folderToDelocalizedUriAndListing.get(path)
            case _          => None
          }) match {
          case x @ Some(_)      => x
          case None if optional => None
          case _                => throw new Exception(s"Did not localize folder ${uri}")
        }
      }
      (t, v) match {
        case (_, f: VFile) =>
          translateFileUri(f.uri).map { newUri =>
            val newSecondaryFiles =
              f.secondaryFiles.flatMap(pathTranslator(_, None, optional = false))
            f.copy(uri = newUri, secondaryFiles = newSecondaryFiles)
          }
        case (Some(TFile), VString(uri))      => translateFileUri(uri).map(VFile(_))
        case (Some(TDirectory), VString(uri)) => translateFolderUri(uri).map(VFolder(_))
        case (_, f: VFolder) =>
          translateFolderUri(f.uri).map(newUri => f.copy(uri = newUri))
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
      case (t, v, _, _) => (t, Value.transform(v, Some(t), pathTranslator))
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
