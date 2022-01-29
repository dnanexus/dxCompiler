package dx.executor

import java.nio.file.{Files, Path}
import dx.api.{DxFile, DxJob, DxPath, FileUpload, InstanceTypeRequest}
import dx.core.getVersion
import dx.core.io.{
  DxdaManifest,
  DxdaManifestBuilder,
  DxfuseManifest,
  DxfuseManifestBuilder,
  StreamFiles
}
import dx.core.ir.{DxName, Type, TypeSerde, Value, ValueSerde}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.{
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
import dx.util.protocols.{DxFileSource, DxFolderSource}
import spray.json._

object TaskExecutor {
  val MaxDisambiguationDirs: Int = 5000
  // functions we can call in the command script
  val DownloadDxda = "download_dxda"
  val DownloadDxfuse = "download_dxfuse"
  val BeforeCommand = "before_command"
  val RunCommand = "run_command"
  val UploadFiles = "upload_files"
  val Cleanup = "cleanup"
}

object TaskExecutorResult extends Enum {
  type TaskExecutorResult = Value
  val Success, Relaunch = Value
}

abstract class TaskExecutor(jobMeta: JobMeta,
                            streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
                            checkInstanceType: Boolean,
                            traceLengthLimit: Int = 10000) {

  private val workerPaths = jobMeta.workerPaths
  private val fileResolver = jobMeta.fileResolver
  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger

  protected def trace(msg: String, minLevel: Int = TraceLevel.Verbose): Unit = {
    logger.traceLimited(msg, traceLengthLimit, minLevel)
  }

  def executorName: String

  private def printDirTree(): Unit = {
    if (logger.traceLevel >= TraceLevel.VVerbose) {
      trace("Directory structure:", TraceLevel.VVerbose)
      val (_, stdout, _) = SysUtils.execCommand("ls -lR", None)
      trace(stdout + "\n", TraceLevel.VVerbose)
    }
  }

  /**
    * Returns the IR type and value for all task input variables, including default values for any
    * missing optional parameters.
    */
  protected def getInputVariables: Map[DxName, (Type, Value)]

  /**
    * Returns the minimal (i.e. cheapest) instance type that is
    * sufficient to the task's resource requirements.
    */
  protected def getInstanceTypeRequest(inputs: Map[DxName, (Type, Value)]): InstanceTypeRequest

  /**
    * Returns a mapping of output field names IR types.
    */
  protected def outputTypes: Map[DxName, Type]

  /**
    * Should we try to stream the file(s) associated with the given input parameter?
    */
  protected def streamFileForInput(parameterName: DxName): Boolean

  /**
    * Returns all schema types referenced in the inputs.
    */
  protected def getSchemas: Map[String, TSchema]

  private class PathsToLocalize {
    var virtualFiles: Vector[StringFileNode] = Vector.empty
    var localFiles: Set[LocalFileSource] = Set.empty
    var filesToStream: Set[AddressableFileSource] = Set.empty
    var filesToDownload: Set[AddressableFileSource] = Set.empty

    private class Handler(stream: Boolean) extends TransformHandler {
      // Updates all `VFolder` values (including those nested within `value`) to
      // set their `listing`s if missing. Also extracts all file `VFile` values into
      // a `PathsToLocalize` object as follows:
      // * If `contents` is non-empty, the file is added to `virtualFiles`
      // * If `address` resolves to a `LocalFileSource`, the file is added to `localFiles`
      // * Otherwise, if `stream` is true, the file is added to `filesToStream`, otherwise
      //   it is added to `filesToDownload`
      def apply(value: Value, t: Option[Type] = None, optional: Boolean = false): Option[Value] = {
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
            apply(path) match {
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

        // gets the listing for a folder - if the folder does not exist, returns an empty listing
        def getFolderListing(uri: String): Vector[FileSource] = {
          val fs =
            try {
              fileResolver.resolveDirectory(uri)
            } catch {
              case _: Throwable => return Vector()
            }
          if (!fs.isDirectory) {
            throw new Exception(s"Directory-type input is not a directory: ${fs}")
          }
          if (fs.exists && fs.isListable) {
            fs.listing(recursive = true)
          } else {
            Vector()
          }
        }

        (t, value) match {
          case (Some(TFile) | None, file: VFile) if file.contents.nonEmpty =>
            virtualFiles :+= StringFileNode(file.contents.get, file.uri)
            Some(file.copy(contents = None, secondaryFiles = updateListing(file.secondaryFiles)))
          case (Some(TFile) | None, file: VFile) =>
            handleFileUri(file.uri)
            Some(file.copy(secondaryFiles = updateListing(file.secondaryFiles)))
          case (Some(TFile), VString(uri)) =>
            handleFileUri(uri)
            Some(VFile(uri))
          case (Some(TDirectory) | None, folder: VFolder) if folder.listing.isEmpty =>
            // fill in the listing if it is not provided
            Some(folder.copy(listing = Some(createListing(getFolderListing(folder.uri)))))
          case (Some(TDirectory) | None, folder: VFolder) =>
            Some(folder.copy(listing = Some(updateListing(folder.listing.get))))
          case (Some(TDirectory), VString(uri)) =>
            val listing = Some(getFolderListing(uri)).filter(_.nonEmpty).map(createListing)
            Some(VFolder(uri, listing = listing))
          case (Some(TDirectory) | None, listing: VListing) =>
            Some(listing.copy(items = updateListing(listing.items)))
          case _ => None
        }
      }
    }

    def updateListingsAndExtractFiles(
        inputs: Map[DxName, (Type, Value)]
    ): Map[DxName, (Type, Value)] = {
      inputs.map {
        case (dxName, (t, v)) =>
          // whether to stream all the files/directories associated with this input
          val stream = streamFiles == StreamFiles.All ||
            (streamFiles == StreamFiles.PerFile && streamFileForInput(dxName))
          val handler = new Handler(stream)
          dxName -> (t, Value.transform(v, Some(t), handler))
      }
    }

    def resolveStaticDependencies(uris: Set[String]): Map[String, (Type, Value)] = {
      val handler = new Handler(stream = false)
      uris.map {
        case uri if uri.endsWith("/") =>
          val resolved = handler.apply(VFolder(uri), Some(TDirectory)) match {
            case Some(f: VFolder) => f
            case other            => throw new Exception(s"expected VFolder, not ${other}")
          }
          uri -> (TDirectory, resolved)
        case uri =>
          // update the URI to match the uri of a raw input to avoid downloading the same file twice
          // e.g. the file is stored in the static dependencies and is also the default value in a
          // field and is not overridden - we'd get the same file twice, but in the first case the
          // uri is `dx://file-xxx::/path/to/file` and in the second it might be
          // `dx://project-xxx:/path/to/file` or `dx://project-xxx:file-yyy`.
          // TODO: a better way to handle duplicate files with different URIs would be to add
          //  another layer of indirection, such that all URIs that resolve to the same file ID
          //  map to the same "cannonical" AddressableFileSource
          val uriToResolve = fileResolver.resolve(uri) match {
            case fs: DxFileSource =>
              DxFile.format(fs.dxFile.id,
                            fs.dxFile.describe().folder,
                            fs.dxFile.describe().name,
                            None)
            case _ => uri
          }
          val resolved = handler.apply(VFile(uriToResolve), Some(TFile)) match {
            case Some(f: VFile) => f
            case other          => throw new Exception(s"expected VFile, not ${other}")
          }
          uri -> (TFile, resolved)
      }.toMap
    }
  }

  private val streamingDir = workerPaths.getDxfuseMountDir().asJavaPath
  private val localizationDirs = Vector(workerPaths.getInputFilesDir().asJavaPath, streamingDir)

  private class InputFinalizer(
      uriToPath: Map[String, Path],
      localizer: SafeLocalizationDisambiguator
  ) extends TransformHandler {
    private var pathToUri: Map[Path, String] = Map.empty
    // mappings to track final paths for files/folders - we want to ensure that
    // identical source VFiles/VFolders map to the same localized path
    private var sourceToFinalFile: Map[VFile, Path] = Map.empty
    private var sourceToFinalFolder: Map[VFolder, Path] = Map.empty

    def finalizeListing(listing: Vector[PathValue], listingParent: Path): Vector[PathValue] = {
      listing.map { value =>
        finalizeValue(value, None, Some(listingParent)) match {
          case Some(p: PathValue) => p
          case other              => throw new Exception(s"expected PathValue, not ${other}")
        }
      }
    }

    def finalizeValue(value: Value, t: Option[Type], parent: Option[Path]): Option[Value] = {
      (t, value) match {
        case (Some(TFile) | None, f: VFile) =>
          val path = uriToPath(f.uri)
          val parentDir = parent.getOrElse(path.getParent)
          lazy val fs = fileResolver.resolve(f.uri)
          val newPath = sourceToFinalFile.get(f) match {
            case Some(newPath) => newPath
            case None
                if f.basename.nonEmpty ||
                  !localizationDirs.exists(path.startsWith) ||
                  (parent.nonEmpty && !path.startsWith(parent.get)) =>
              // Either we are localizing into a specific directory, the path is a local or a
              // virtual file, or it is a remote file that needs its name changed - use the
              // localizer to determine the new path and then create a symbolic link.
              val pathToLocalize = f.basename.map(parentDir.resolve).getOrElse(path)
              val newPath = localizer
                .getLocalPath(fileResolver.fromPath(pathToLocalize, isDirectory = Some(false)),
                              parent)
              Files.createLink(newPath, path)
              sourceToFinalFile += (f -> newPath)
              newPath
            case None if path.startsWith(streamingDir) && localizer.getTargetDir(fs).isDefined =>
              // we have streaming and non-streaming files from the same source container -
              // link the streaming file into the download dir
              val newPath = localizer.getLocalPath(fs)
              Files.createLink(newPath, path)
              sourceToFinalFile += (f -> newPath)
              newPath
            case _ => path
          }
          pathToUri += (newPath -> f.uri)
          // recursively resolve the secondary files, linking them to be adjacent to the main file
          val newSecondaryFiles = finalizeListing(f.secondaryFiles, newPath.getParent)
          Some(f.copy(uri = newPath.toString, secondaryFiles = newSecondaryFiles))
        case (Some(TDirectory) | None, f: VFolder) if f.listing.isDefined =>
          // link all the files/subdirectories in the listing into the parent folder
          val dirPath = sourceToFinalFolder.getOrElse(
              f, {
                val fs = (fileResolver.resolveDirectory(f.uri), f.basename) match {
                  case (fs, Some(basename)) =>
                    fs.getParent
                      .map(_.resolveDirectory(basename))
                      .getOrElse(throw new Exception("cannot rename root directory"))
                  case (fs, None) => fs
                }
                val dirPath = localizer.getLocalPath(fs, parent)
                FileUtils.createDirectories(dirPath)
                sourceToFinalFolder += (f -> dirPath)
                dirPath
              }
          )
          pathToUri += (dirPath -> f.uri)
          val newListing = finalizeListing(f.listing.get, dirPath)
          Some(f.copy(uri = dirPath.toString, listing = Some(newListing)))
        case (Some(TDirectory) | None, l: VListing) =>
          // a listing is a virtual directory and is unique (i.e. we never want two
          // listings pointing at the same physical directory), so we source it from
          // a temp dir
          val tempDir = Files.createTempDirectory("listing")
          tempDir.toFile.deleteOnExit()
          val fs = fileResolver.fromPath(tempDir.resolve(l.basename), isDirectory = Some(true))
          val dirPath = localizer.getLocalPath(fs, parent)
          FileUtils.createDirectories(dirPath)
          val newListing = finalizeListing(l.items, dirPath)
          Some(l.copy(items = newListing))
        case _ => None
      }
    }

    def apply(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      finalizeValue(value, t, None)
    }

    def finalizeInputs(
        inputs: Map[DxName, (Type, Value)]
    ): (Map[DxName, (Type, Value)], Map[Path, String]) = {
      val finalizedInputs = inputs.map {
        case (dxName, (t, v)) =>
          dxName -> (t, Value.transform(v, Some(t), this))
      }
      (finalizedInputs, pathToUri)
    }

    def finalizeStaticDependencies(
        deps: Map[String, (Type, Value)]
    ): Map[String, (Type, Value)] = {
      deps.map {
        case (uri, (t, v)) =>
          uri -> (t, Value.transform(v, Some(t), this))
      }
    }
  }

  /**
    * Generates and writes command script(s) to disk.
    * @param localizedInputs task inputs with localized files
    * @param localizedDependencies mapping of URI to (type, value) for any static
    *                              dependencies that may need to be updated in
    *                              the source code.
    * @return (updatedInputs, hasCommand, successCodes) where updated
    *         inputs is localizedInputs updated with any additional (non-input)
    *         variables that may be required to evaluate the outputs, hasCommand
    *         is whether a command script was actually written, and successCodes
    *         are the return codes that indicate the script was executed
    *         successfully. If successCodes is None, then all return codes are
    *         considered successful.
    */
  protected def writeCommandScript(
      localizedInputs: Map[DxName, (Type, Value)],
      localizedDependencies: Option[Map[String, (Type, Value)]]
  ): (Boolean, Option[Set[Int]])

  /**
    * Evaluates the outputs of the task. Returns mapping of output parameter
    * name to (type, value) and to (Set of tags, and Map of properties), where tags
    * and properites only apply to output files (or collections thereof).
    * @param localizedInputs the job inputs, with files localized to the worker
    */
  protected def evaluateOutputs(
      localizedInputs: Map[DxName, (Type, Value)]
  ): (Map[DxName, (Type, Value)], Map[DxName, (Set[String], Map[String, String])])

  private class Listings(localizedInputs: Map[DxName, (Type, Value)]) {
    // lazily build mapping of input path to VFolder
    private object ListingWalkHandler extends WalkHandler[Map[Path, VFolder]] {
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

      def apply(value: Value,
                t: Option[Type],
                optional: Boolean,
                folders: Map[Path, VFolder]): Option[Map[Path, VFolder]] = {
        (t, value) match {
          case (Some(TDirectory) | None, f: VFolder) => Some(handlePath(f, folders))
          case _                                     => None
        }
      }
    }

    lazy val localPathToInputFolder: Map[Path, VFolder] = {
      localizedInputs.foldLeft(Map.empty[Path, VFolder]) {
        case (accu, (_, (t, v))) => Value.walk(v, Some(t), accu, ListingWalkHandler)
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
        listing.flatMap {
          case f: VFile =>
            fileResolver.resolve(f.uri) match {
              case local: LocalFileSource =>
                Set(local.canonicalPath) ++ getListingPaths(f.secondaryFiles)
              case _ => Set.empty[Path]
            }
          case f: VFolder =>
            fileResolver.resolveDirectory(f.uri) match {
              case local: LocalFileSource => Set(local.canonicalPath)
              case _                      => Set.empty[Path]
            }
          case _ => Set.empty[Path]
        }
      }.toSet

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
              case other => throw new Exception(s"unexpected item ${other}")
            }
          }

          outputListing.forall(itemEqual)
      }
    }

    // creates a listing of VFile/VFolders from a nested listing of FileSources
    private def createListing(value: Vector[FileSource]): Vector[PathValue] = {
      value.map {
        case fs: LocalFileSource if fs.isDirectory =>
          VFolder(fs.address, listing = Some(createListing(fs.listing())))
        case fs: LocalFileSource => VFile(fs.address)
        case other               => throw new Exception(s"unexpected output listing item ${other}")
      }
    }

    // update outputs - make sure all folders have listings
    private object ListingTransformHandler extends TransformHandler {
      def apply(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
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
    }

    def finalizeOutputs(outputs: Map[DxName, (Type, Value)]): Map[DxName, (Type, Value)] = {
      outputs.map {
        case (dxName, (t, v)) => dxName -> (t, Value.transform(v, Some(t), ListingTransformHandler))
      }
    }
  }

  private class Uploads(tags: Set[String], properties: Map[String, String]) {
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
      val localPath = workerPaths.getVirtualFilesDir(ensureExists = true).resolve(uri).asJavaPath
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

  private class FileExtractor(
      finalizedOutputs: Map[DxName, (Type, Value)],
      tagsAndProperties: Map[DxName, (Set[String], Map[String, String])],
      uriToSourcePath: Map[String, Path],
      localPathToUri: Map[Path, String],
      listings: Listings
  ) extends WalkHandler[Uploads] {
    private val defaultDestParent = if (jobMeta.useManifests) {
      DxFolderSource.ensureEndsWithSlash(jobMeta.manifestProjectAndFolder)
    } else {
      "/"
    }

    // does a local path correspond to a platform file?
    def isDxPath(path: Path): Boolean = {
      localPathToUri.get(path).exists { uri =>
        val sourcePath = uriToSourcePath(uri)
        localizationDirs.exists(sourcePath.startsWith)
      }
    }

    // extract all files to upload
    def apply(value: Value, t: Option[Type], optional: Boolean, ctx: Uploads): Option[Uploads] = {
      def handlePath(innerValue: PathValue,
                     innerType: Option[Type],
                     parent: Option[String] = None): Unit = {
        (innerType, innerValue) match {
          case (Some(TFile) | None, f: VFile) if f.contents.isDefined =>
            // a file literal
            val name = f.basename.getOrElse(f.uri)
            val destFolder = parent.getOrElse(defaultDestParent)
            ctx.addVirtualFile(f.contents.get, f.uri, DxFolderSource.join(destFolder, name))
            f.secondaryFiles.foreach { path =>
              handlePath(path, None, Some(destFolder))
            }
          case (Some(TFile) | None, f: VFile) =>
            fileResolver.resolve(f.uri) match {
              case fs if optional && !fs.exists => ()
              case fs if !fs.exists =>
                throw new Exception(s"non-optional output file ${fs} does not exist")
              case local: LocalFileSource
                  if parent.isDefined || f.basename.isDefined || !isDxPath(local.canonicalPath) =>
                // Either the file was created by the command script, it originated as
                // a local or virtual file, it needs to be uploaded with a new name, or
                // we are uploading a listing/directory and so need to make a new copy
                // of the file regardless of its origin
                val name = f.basename.getOrElse(local.canonicalPath.getFileName.toString)
                val destFolder = parent.getOrElse(defaultDestParent)
                ctx.addFile(local.canonicalPath, DxFolderSource.join(destFolder, name))
                f.secondaryFiles.foreach { path =>
                  handlePath(path, None, Some(destFolder))
                }
              case _ => ()
            }
          case (Some(TDirectory) | None, f: VFolder) =>
            fileResolver.resolveDirectory(f.uri) match {
              case fs if optional && !fs.exists => ()
              case fs if !fs.exists =>
                throw new Exception(s"non-optional output directory ${fs} does not exist")
              case local: LocalFileSource
                  if parent.isDefined ||
                    f.basename.isDefined ||
                    !isDxPath(local.canonicalPath) ||
                    f.listing.exists { listing =>
                      !listings.listingsEqual(local.canonicalPath, listing)
                    } =>
                // Either the directory was created by the command script, it originated as
                // a local directory, it needs to be uploaded with a new name, its listing
                // is different, or we are uploading a subdirectory of a listing/directory
                // and so need to make a new copy of the this directory regardless of its origin
                val name = f.basename.getOrElse(local.canonicalPath.getFileName.toString)
                val destFolder =
                  DxFolderSource.join(parent.getOrElse(defaultDestParent), name, isFolder = true)
                ctx.addDirectory(local.canonicalPath, destFolder)
                f.listing.get.foreach { item =>
                  handlePath(item, None, Some(destFolder))
                }
              case _ => ()
            }
          case (Some(TDirectory) | None, l: VListing) =>
            val destFolder =
              DxFolderSource.join(parent.getOrElse(defaultDestParent), l.basename, isFolder = true)
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

    def extractFiles: (Iterable[FileUpload], Map[String, Path], Map[Path, String]) = {
      val (filesToUpload, virtualFiles, directories) = finalizedOutputs.map {
        case (name, (t, v)) =>
          val (tags, properties) =
            tagsAndProperties.getOrElse(name, (Set.empty[String], Map.empty[String, String]))
          Value.walk(v, Some(t), new Uploads(tags, properties), this).asTuple
      }.unzip3
      (filesToUpload.flatten, virtualFiles.flatten.toMap, directories.flatten.toMap)
    }
  }

  private def delocalizePaths(outputs: Map[DxName, (Type, Value)],
                              uploadedFiles: Map[Path, DxFile],
                              uploadedVirtualFiles: Map[String, Path],
                              uploadedDirectories: Map[Path, String],
                              localPathToUri: Map[Path, String],
                              listings: Listings): Map[DxName, (Type, Value)] = {
    // For folders/listings, we have to determine the project output folder, because the platform
    // won't automatically update directory URIs like it does for data object links. The output
    // folder is either the manifest folder or the concatination of the job/analysis output folder
    // and the container output folder.
    val parentProjectFolder = DxFolderSource.ensureEndsWithSlash(jobMeta.projectOutputFolder)

    // replace local paths with remote URIs
    def handlePath(value: PathValue,
                   t: Option[Type] = None,
                   optional: Boolean = false,
                   parent: Option[String] = None): Value = {
      (t, value) match {
        case (Some(TFile) | None, f: VFile) if f.contents.isDefined =>
          // a virtual file that was uploaded - replace the URI with the dx URI and
          // unset the contents
          val dxFile = uploadedFiles(uploadedVirtualFiles(f.uri))
          VFile(dxFile.asUri, secondaryFiles = f.secondaryFiles.map(handlePath(_)).map {
            case p: PathValue => p
            case other        => throw new Exception(s"not a PathValue ${other}")
          })
        case (Some(TFile) | None, f: VFile) =>
          val newSecondaryFiles = f.secondaryFiles.map(handlePath(_)).map {
            case p: PathValue => p
            case other        => throw new Exception(s"not a PathValue ${other}")
          }
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
              val folderParent = parent.getOrElse(parentProjectFolder)
              val name = f.basename.getOrElse(local.canonicalPath.getFileName.toString)
              val folder = DxFolderSource.join(folderParent, name)
              VFolder(
                  DxFolderSource.format(dxApi.currentProject.get, folder),
                  basename = Some(name),
                  listing = f.listing.map(_.map(handlePath(_, parent = Some(folder))).map {
                    case p: PathValue => p
                    case other        => throw new Exception(s"not a PathValue ${other}")
                  })
              )
            case local: LocalFileSource
                if localPathToUri.contains(local.canonicalPath) && parent.isEmpty =>
              // an output folder that was an input folder
              VFolder(localPathToUri(local.canonicalPath),
                      listing = listings.getInputListing(local.canonicalPath))
            case _: LocalFileSource =>
              throw new Exception(s"folder was not delocalized: ${f}")
            case _ => f
          }
        case (Some(TDirectory) | None, l: VListing) =>
          val folderParent = parent.getOrElse(parentProjectFolder)
          val folder = DxFolderSource.join(folderParent, l.basename)
          VFolder(
              DxFolderSource.format(dxApi.currentProject.get, folder),
              basename = Some(l.basename),
              listing = Some(l.items.map(handlePath(_, parent = Some(folder))).map {
                case p: PathValue => p
                case other        => throw new Exception(s"not a PathValue ${other}")
              })
          )
        case _ =>
          throw new Exception(s"Could not delocalize ${value} as ${t}")
      }
    }

    def delocalizePaths(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      (t, value) match {
        case (Some(TFile) | Some(TDirectory) | None, p: PathValue) =>
          Some(handlePath(p, t, optional))
        case _ => None
      }
    }

    outputs.map {
      case (dxName, (t, v)) =>
        dxName -> (t, Value.transform(v, Some(t), delocalizePaths))
    }
  }

  private def logFields(fields: Map[DxName, (Type, Value)], description: String): Unit = {
    if (logger.isVerbose) {
      val inputStr = fields
        .map {
          case (dxName, (t, v)) =>
            s"${dxName}: (${TypeSerde.toString(t)}, ${ValueSerde.toString(v, verbose = true)})"
        }
        .mkString("\n  ")
      trace(s"${description}:\n  ${inputStr}")
    }
  }

  /**
    * Executes the task. This implements the full lifecycle of task execution:
    * 1. If the instance type is determined dynamically, computes the instance
    *    requirements and checks whether the current instance type meets/exceeds
    *    those requirements. If not, relaunches the task with a new instance type.
    * 2. Evaluate inputs.
    * 3. Localize input files.
    * 4. Evaluate and run the command script.
    * 5. Evaluate outputs.
    * 6. Delocalize output files.
    * @return true if the task was executed successfully, or false if was relaunched
    */
  def apply(): TaskExecutorResult.TaskExecutorResult = {
    // setup the utility directories that the task-runner employs
    workerPaths.createCleanDirs()

    if (logger.isVerbose) {
      trace(s"dxCompiler version: ${getVersion}")
      trace(s"debugLevel=${logger.traceLevel}")
      trace(s"Task source code:\n${jobMeta.sourceCode}", traceLengthLimit)
      printDirTree()
    }

    // Evaluates input values and makes sure all VFolder inputs have their `listing`s set.
    // We do this to ensure that only the files present in the folder at the beginning of
    // task execution are localized into the inputs folder, and any changes to the
    // directory during runtime aren't synchronized to the worker. This step also extracts
    // all the files from the inputs that need to be localized. Input files and directories
    // are represented as URIs (dx://proj-xxxx:file-yyyy::/A/B/C.txt).
    // TODO: it would be nice to extract dx:// links from VString values - this will happen
    //  in the case where the container is a dx file and being passed in as an input
    //  parameter - so that they could be downloaded using dxda. However, this would also
    //  require some way for the downloaded image tarball to be discovered and loaded. For now,
    //  we rely on DockerUtils to download the image (via DxFileSource, which uses the dx API
    //  to download the file).
    val pathsToLocalize = new PathsToLocalize
    val inputs = pathsToLocalize.updateListingsAndExtractFiles(getInputVariables)

    if (checkInstanceType) {
      // calculate the required instance type
      trace("Computing instance type to request")
      val instanceTypeRequest: InstanceTypeRequest = getInstanceTypeRequest(inputs)
      trace(s"Requesting instance type: $instanceTypeRequest")
      val requestedInstanceType = jobMeta.instanceTypeDb.apply(instanceTypeRequest).name
      trace(s"Requested instance type: ${requestedInstanceType}")

      // Check if we are already on a sufficent instance type. We relaunch if
      // the current instance type is not sufficient for the task requirements.
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
        // launch a sub-job with the same inputs and the dynamically calculated instance type
        val dxSubJob: DxJob = dxApi.runSubJob(
            "body",
            Some(requestedInstanceType),
            JsObject(jobMeta.rawJsInputs.map {
              case (dxName, jsv) => dxName.encoded -> jsv
            }),
            Vector.empty,
            jobMeta.delayWorkspaceDestruction
        )
        jobMeta.writeExecutionOutputLinks(dxSubJob, outputTypes)
        return TaskExecutorResult.Relaunch
      }
    }

    // If the task has any hard-coded remote paths, localize them and later update the source
    // code to replace the remote paths with local ones
    val staticDependencies = Option.when(jobMeta.fileDependencies.nonEmpty)(
        pathsToLocalize.resolveStaticDependencies(jobMeta.fileDependencies)
    )

    // localize virtual files to /home/dnanexus/virtual
    trace(s"Virtual files = ${pathsToLocalize.virtualFiles}")
    val virtualUriToPath = pathsToLocalize.virtualFiles.map { fs =>
      val localPath =
        fs.localizeToDir(workerPaths.getVirtualFilesDir(ensureExists = true).asJavaPath)
      fs.name -> localPath
    }.toMap

    // local files already have a path
    // TODO: should we use originalPath rather than address?
    val localUriToPath = pathsToLocalize.localFiles.map(fs => fs.address -> fs.canonicalPath).toMap

    // This object handles mapping FileSources to local paths and deals with file name
    // collisions in the manner specified by the WDL spec. All remote input files except those
    // streamed by dxfuse are placed in subfolders of the /home/dnanexus/inputs directory.
    val localizer = SafeLocalizationDisambiguator.create(
        rootDir = workerPaths.getInputFilesDir().asJavaPath,
        separateDirsBySource = true,
        createDirs = true,
        disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
        logger = logger
    )

    // Create manifests for dxfuse and/or dxda to localize all input files. Extracted files are
    // added to on or the other manifest based on the value of `streamFiles` and whether streaming
    // has been enabled individually for any of the inputs.
    val downloadUriToPath = if (pathsToLocalize.filesToDownload.nonEmpty) {
      // Build dxda manifest to localize all non-streaming remote files
      trace(s"Files to download: ${pathsToLocalize.filesToDownload}")
      val downloadFileSourceToPath: Map[AddressableFileSource, Path] =
        localizer.getLocalPaths(pathsToLocalize.filesToDownload)
      val downloadUriToPath = downloadFileSourceToPath.map {
        case (fs, path) => fs.address -> path
      }
      DxdaManifestBuilder(dxApi, logger)
        .apply(downloadFileSourceToPath.collect {
          case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
        })
        .foreach {
          case DxdaManifest(manifestJs) =>
            // write the manifest to a file
            FileUtils.writeFileContent(workerPaths.getDxdaManifestFile().asJavaPath,
                                       manifestJs.prettyPrint)
            // run dxda via a subprocess
            jobMeta.runJobScriptFunction(TaskExecutor.DownloadDxda)
        }
      downloadUriToPath
    } else {
      Map.empty
    }

    val streamUriToPath = if (pathsToLocalize.filesToStream.nonEmpty) {
      // build dxfuse manifest to localize all straming remote files and folders
      trace(s"Files to stream: ${pathsToLocalize.filesToStream}")
      // use a different localizer for the dxfuse mount point
      val streamingLocalizer = SafeLocalizationDisambiguator.create(
          rootDir = workerPaths.getDxfuseMountDir().asJavaPath,
          existingPaths = localizer.getLocalizedPaths,
          separateDirsBySource = true,
          createDirs = false,
          disambiguationDirLimit = TaskExecutor.MaxDisambiguationDirs,
          logger = logger
      )
      val streamFileSourceToPath: Map[AddressableFileSource, Path] =
        streamingLocalizer.getLocalPaths(pathsToLocalize.filesToStream)
      val streamUriToPath = streamFileSourceToPath.map {
        case (fs, path) => fs.address -> path
      }
      DxfuseManifestBuilder(workerPaths, dxApi, logger)
        .apply(streamFileSourceToPath.collect {
          case (dxFs: DxFileSource, localPath) => dxFs.dxFile -> localPath
        })
        .foreach {
          case DxfuseManifest(manifestJs) =>
            FileUtils.writeFileContent(workerPaths.getDxfuseManifestFile().asJavaPath,
                                       manifestJs.prettyPrint)
            // run dxfuse via a subprocess
            jobMeta.runJobScriptFunction(TaskExecutor.DownloadDxfuse)
        }
      streamUriToPath
    } else {
      Map.empty
    }

    val uriToSourcePath =
      (virtualUriToPath ++ localUriToPath ++ downloadUriToPath ++ streamUriToPath).flatMap {
        case (uri, path) =>
          // add alternate forms of dx URIs that may be needed for lookup later
          fileResolver.resolve(uri) match {
            case fs: DxFileSource =>
              Set(
                  Some(uri),
                  Some(s"${DxPath.DxUriPrefix}${fs.dxFile.id}"),
                  fs.dxFile.project.map(proj => s"${DxPath.DxUriPrefix}${proj.id}:${fs.dxFile.id}")
              ).flatten.map(uri => uri -> path)
            case _ => Set(uri -> path)
          }
      }

    // Finalize all input files, folders, and listings.
    // - For a file that is not in the context of a folder:
    //   - If it is a local file, we link it into a disambiguation dir, naming it with its basename
    //     if it has one.
    //   - Otherwise, if it has a basename, we create a link with the new name in the same folder to
    //     the original file, throwing an exception if there is a naming collision.
    //   - If it has secondary files, they are linked into the same directory as the main file,
    //     creating any subfolders, and throwing an exception if there is a naming collision.
    // - For a directory or listing, we create a new disambiguation dir and recursively link in all
    //   the files it contains, creating any subfolders.
    // - For streaming files, if the source container is the same as one managed by download
    //   localizer, link the streaming file into the download directory.
    val inputFinalizer = new InputFinalizer(uriToSourcePath, localizer)
    val (localizedInputs, localPathToUri) = inputFinalizer.finalizeInputs(inputs)
    logFields(localizedInputs, "Localized inputs")

    // Finalize any static dependencies and update the task source code if necessary
    val finalizedDependencies = staticDependencies.map(inputFinalizer.finalizeStaticDependencies)

    // make the inputs folder read-only
    jobMeta.runJobScriptFunction(TaskExecutor.BeforeCommand)

    // Evaluate the command script and writes it to disk. Inputs are supplemented with any local
    // file paths created when evaluating the command script and are serialized for use in the next
    // phase.
    val (hasCommand, successCodes) = writeCommandScript(localizedInputs, finalizedDependencies)
    if (hasCommand) {
      // run the command script
      jobMeta.runJobScriptFunction(TaskExecutor.RunCommand, successCodes, forwardStd = true)
    }

    // evaluate output expressions
    val (localizedOutputs, tagsAndProperties) = evaluateOutputs(localizedInputs)
    logFields(localizedOutputs, "Localized outputs")

    // Upload output files/directories and replace local paths with remote URIs in the output
    // values. An output file/folder may have been created on the worker or it may be an input that
    // originated either on the worker or from the platform. We don't want to reupload
    // files/directories that already have an identical copy at the specified output location on
    // the platform, but otherwise we need to upload them.
    if (localizedOutputs.nonEmpty) {
      if (logger.isVerbose) {
        trace(s"Delocalizing outputs")
        printDirTree()
      }

      val listings = new Listings(localizedInputs)

      // set listings on output directories
      val finalizedOutputs = listings.finalizeOutputs(localizedOutputs)

      // get files to upload
      val fileExtractor = new FileExtractor(finalizedOutputs,
                                            tagsAndProperties,
                                            uriToSourcePath,
                                            localPathToUri,
                                            listings)
      val (filesToUpload, virtualFiles, directories) = fileExtractor.extractFiles

      // upload files
      // TODO: handle output parameter tags and properties
      val uploadedFiles = if (filesToUpload.nonEmpty) {
        if (logger.isVerbose) {
          trace(s"Uploading files:\n  ${filesToUpload.map(_.source).mkString("\n")}")
        }
        // upload all files in parallel
        jobMeta.uploadFiles(filesToUpload)
        // TODO: use dxua or dxfuse instead
        // jobMeta.runJobScriptFunction(TaskExecutor.UploadFiles)
      } else {
        Map.empty[Path, DxFile]
      }

      // replace local paths with remote paths in outputs
      val delocalizedOutputs = delocalizePaths(finalizedOutputs,
                                               uploadedFiles,
                                               virtualFiles,
                                               directories,
                                               localPathToUri,
                                               listings)

      logFields(delocalizedOutputs, "Delocalized outputs")

      // serialize the outputs to the job output file
      jobMeta.writeOutputs(delocalizedOutputs)
    }

    // peform any cleanup - unmount dxfuse dir
    jobMeta.runJobScriptFunction(TaskExecutor.Cleanup)

    TaskExecutorResult.Success
  }
}
