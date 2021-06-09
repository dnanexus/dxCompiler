package dx.core.ir

import java.nio.charset.Charset
import java.nio.file.{Files, Path, Paths}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.{FileUtils, JsUtils, Logger, SysUtils}
import spray.json._

import scala.collection.immutable.TreeSeqMap

/**
  * Interface for mounting or appending to a squashfs image.
  * @param image path to the image file
  * @param mountDir directory where to mount the file system - if Right, is
  *                 used as the actual mount point; if Left, is the parent
  *                 directory of a randomly named mount point; if null, a
  *                 random mount point is created in tmp
  * @param logger Logger
  */
case class SquashFs(image: Path)(
    mountDir: Option[Either[Path, Path]] = None,
    alreadyMounted: Boolean = false,
    logger: Logger = Logger.get
) {
  lazy val mountPoint: Path = {
    mountDir match {
      case Some(Right(dir)) => dir
      case Some(Left(dir)) =>
        Files.createTempDirectory(dir, image.getFileName.toString)
      case None =>
        Files.createTempDirectory(image.getFileName.toString)
    }
  }
  private var mounted: Boolean = alreadyMounted
  // whether the archive file is a placeholder that needs to be deleted
  // before writing the archive
  private var newFile: Boolean = Files.exists(image) && Files.size(image) == 0
  private lazy val isLinux = System.getProperty("os.name").toLowerCase.startsWith("linux")

  def isMounted: Boolean = mounted

  private def mount(sudo: Boolean): Unit = {
    val prefix = if (sudo) "sudo -n " else ""
    if (isLinux) {
      SysUtils.execCommand(
          s"${prefix}mount -t squashfs -o loop ${image.toString} ${mountPoint.toString}"
      )
    } else {
      SysUtils.execCommand(s"${prefix}squashfuse ${image.toString} ${mountPoint.toString}")
    }
    mounted = true
  }

  /**
    * Mounts the archive if it has not already been mounted.
    */
  def mount(): Unit = {
    if (!mounted) {
      try {
        // first try without sudo
        mount(sudo = false)
      } catch {
        case _: Throwable =>
          try {
            // try again with sudo
            mount(sudo = true)
          } catch {
            case ex: Throwable =>
              throw new RuntimeException(s"unable to mount archive ${image} at ${mountPoint}", ex)
          }
      }
    }
  }

  private def unmount(sudo: Boolean): Unit = {
    val prefix = if (sudo) "sudo -n " else ""
    val forceOpt = if (isLinux) "-l" else "-f"
    SysUtils.execCommand(s"${prefix}umount ${forceOpt} ${mountPoint.toString}")
    mounted = false
  }

  /**
    * Unmounts the archive if it has previously been mounted.
    */
  def unmount(): Unit = {
    if (isMounted) {
      try {
        // first try without sudo
        unmount(sudo = false)
      } catch {
        case _: Throwable =>
          // try again with sudo
          try {
            unmount(sudo = true)
          } catch {
            case ex: Throwable =>
              logger.error(s"error unmounting archive at ${mountPoint}", Some(ex))
          }
      }
    }
  }

  /**
    * Resolves the given relative path within the mounted directory.
    */
  def resolve(path: String): Path = {
    mountPoint.resolve(path)
  }

  /**
    * Resolves the given relative path within the mounted directory.
    */
  def resolve(path: Path): Path = {
    mountPoint.resolve(path)
  }

  private lazy val appendDir: Path = Files.createTempDirectory("append")

  /**
    * Appends a file to the archive. Throws an exception if the archive has already been
    * mounted. We *could* simply unmount and remount the file system after the append,
    * but this might cause problems for another thread/process trying to access the file
    * system.
    * @param files map from absolute source path to target name
    * @param subdir optional subdir to place the files in within the archive
    * @param removeSource remove the original file after it has been added to the
    *                       archive
    */
  def appendAll(files: Vector[Path],
                subdir: Option[Path] = None,
                removeSource: Boolean = false): Unit = {
    if (isMounted) {
      throw new RuntimeException("cannot append to an archive that is mounted")
    }
    val opts = if (newFile) {
      // overwrite the new file
      newFile = false
      "-noappend"
    } else {
      ""
    }
    try {
      if (subdir.isDefined) {
        // mksquashfs can't easily add a specific file in a subdir, so we need to
        // create a temp directory with the desired sub-directory structure, then
        // create hard links to the files within the leaf dir, then append the files,
        // and finally delete the links
        val lnDir = appendDir.resolve(subdir.get)
        SysUtils.execCommand(
            s"""mkdir -p ${lnDir}
               |&& ln ${files.mkString(" ")} ${lnDir}
               |&& mksquashfs ${appendDir.toString} ${image.toString} ${opts}
               |&& rm -Rf ${appendDir}/*""".stripMargin.replaceAll("\n", " ")
        )
      } else {
        SysUtils.execCommand(
            s"mksquashfs ${files.mkString(" ")} ${image.toString} ${opts}"
        )
      }
    } catch {
      case ex: Throwable =>
        throw new RuntimeException(
            s"error adding files (${files.mkString(",")}) to archive ${image}",
            ex
        )
    }
    if (removeSource) {
      files.foreach {
        case srcPath if Files.isWritable(srcPath) =>
          try {
            Files.delete(srcPath)
          } catch {
            case ex: Throwable =>
              logger.error(s"error deleting ${srcPath} after adding it to archive ${image}",
                           Some(ex))
          }
      }
    }
  }

  def append(file: Path, subdir: Option[Path] = None, removeSource: Boolean = false): Unit = {
    appendAll(Vector(file), subdir, removeSource)
  }

  /**
    * Note that it is only safe to use this for temporary serialization within
    * the same worker session.
    */
  def toJson: JsValue = {
    val (mountPointJs, mountDirJs) = if (mounted) {
      (JsString(mountPoint.toString), JsNull)
    } else {
      mountDir match {
        case Some(Left(path))  => (JsNull, JsString(path.toString))
        case Some(Right(path)) => (JsString(path.toString), JsNull)
        case _                 => (JsNull, JsNull)
      }
    }
    JsObject(
        "image" -> JsString(image.toString),
        "mountPoint" -> mountPointJs,
        "mountDir" -> mountDirJs,
        "mounted" -> JsBoolean(mounted)
    )
  }
}

object SquashFs {
  def fromJson(jsValue: JsValue): SquashFs = {
    jsValue.asJsObject.getFields("image", "mountPoint", "mountDir", "mounted") match {
      case Seq(imageJs, mountPointJs, mountDirJs, mountedJs) =>
        val image = Paths.get(JsUtils.getString(imageJs))
        val mounted = JsUtils.getBoolean(mountedJs)
        val mountDir = if (mounted) {
          Some(Right(Paths.get(JsUtils.getString(mountPointJs))))
        } else {
          (mountDirJs, mountPointJs) match {
            case (JsNull, JsString(path)) => Some(Right(Paths.get(path)))
            case (JsString(path), JsNull) => Some(Left(Paths.get(path)))
            case (JsNull, JsNull)         => None
            case other =>
              throw new Exception(s"invalid mountPoint and mountDir values ${other}")
          }
        }
        SquashFs(image)(mountDir, mounted)
      case _ =>
        throw new Exception(s"invalid serialized SquashFs value ${jsValue}")
    }
  }
}

/**
  * An Archive is a squashfs image that contains 1) a JSON file (called manifest.json)
  * with the serialized representation of a complex value, along with its serialized
  * type information, and 2) the files referenced by all of the VFile values nested
  * within the complex value.
  * TODO: add option to support compressing the file system
  */
trait Archive {
  val path: Path
  val irType: Type
  val irValue: Value
  val localized: Boolean
}

object Archive {
  val ArchiveFilePrefix = "archive___"
  val ArchiveFileSuffix = ".img"
  val ManifestFile: String = "manifest.json"
  val ManifestValueKey: String = "value"

  def isArchiveFile(name: String): Boolean = {
    name.startsWith(ArchiveFilePrefix) && name.endsWith(ArchiveFileSuffix)
  }

  /**
    * Transforms the paths of File-typed values, which may be contained in
    * a (possibly nested) collection.
    * @param irValue the WdlValue to transform
    * @param irType the WdlType of the value
    * @param transformer function to transform one Path to another
    * @return (updatedValue, pathMap), where the updatedValue is identical
    *         to `irValue` except with all paths of file-typed values updated,
    *         and pathMap is a mapping from old to new paths
    */
  def transformPaths(irValue: Value,
                     irType: Type,
                     transformer: Path => Path): (Value, Map[Path, Path]) = {
    def transformFile(path: String): (VFile, Map[Path, Path]) = {
      val oldPath = Paths.get(path)
      val newPath = transformer(oldPath)
      (VFile(newPath.toString), TreeSeqMap(oldPath -> newPath))
    }
    def transformNoType(innerValue: Value): (Value, Map[Path, Path]) = {
      innerValue match {
        case VFile(s) => transformFile(s)
        case VArray(items) =>
          val (transformedItems, paths) = items.map(transformNoType).unzip
          (VArray(transformedItems), paths.flatten.toMap)
        case VHash(fields) =>
          val (transformedFields, paths) = fields.map {
            case (k, v) =>
              val (value, paths) = transformNoType(v)
              (k -> value, paths)
          }.unzip
          (VHash(transformedFields.to(TreeSeqMap)), paths.flatten.toMap)
        case _ =>
          (innerValue, Map.empty)
      }
    }
    def transformWithType(innerValue: Value, innerType: Type): (Value, Map[Path, Path]) = {
      (innerType, innerValue) match {
        case (TFile, VFile(s))   => transformFile(s)
        case (TFile, VString(s)) => transformFile(s)
        case (_: TCollection, _: VFile) =>
          throw new RuntimeException("nested archive values are not allowed")
        case (TOptional(t), value) if value != VNull =>
          val (v, paths) = transformWithType(value, t)
          (v, paths)
        case (TArray(itemType, _), VArray(items)) =>
          val (transformedItems, paths) = items.map(transformWithType(_, itemType)).unzip
          (VArray(transformedItems), paths.flatten.toMap)
        case (TSchema(name, fieldTypes), VHash(fields)) =>
          val (transformedFields, paths) = fieldTypes.collect {
            case (name, t) if fields.contains(name) =>
              val (transformedValue, paths) = transformWithType(fields(name), t)
              (name -> transformedValue, paths)
          }.unzip
          (VHash(transformedFields.to(TreeSeqMap)), paths.flatten.toMap)
        case (THash, _: VHash) =>
          transformNoType(innerValue)
        case _ =>
          (innerValue, Map.empty)
      }
    }
    transformWithType(irValue, irType)
  }
}

/**
  * An Archive that is not localized.
  * @param path path to the archive (.img) file
  * @param encoding character encoding
  * @param typeAliases type aliases
  * @param packedTypeAndValue the already localized type and value
  * @param mountDir parent dir in which to create the randomly named mount point
  */
case class PackedArchive(path: Path, encoding: Charset = FileUtils.DefaultEncoding)(
    typeAliases: Map[String, TSchema] = Map.empty,
    packedTypeAndValue: Option[(Type, Value)] = None,
    mountDir: Option[Path] = None
) extends Archive {
  if (!Files.exists(path)) {
    throw new Exception(s"${path} does not exist")
  } else if (Files.isDirectory(path)) {
    throw new Exception(s"${path} is not a file")
  }

  override val localized: Boolean = false
  private val archive: SquashFs = SquashFs(path)(mountDir.map(Left(_)))
  private lazy val manifestPath: Path = archive.resolve(Archive.ManifestFile)

  private def readManifest: (Type, Value) = {
    if (!archive.isMounted) {
      archive.mount()
      sys.addShutdownHook({
        if (archive.isMounted) {
          archive.unmount()
        }
      })
    }
    val manifestJs =
      try {
        JsUtils.jsFromFile(manifestPath).asJsObject
      } catch {
        case ex: Throwable =>
          throw new RuntimeException(s"unable to read archive manifest ${manifestPath}", ex)
      }
    val (irType, _) = TypeSerde.deserializeOne(manifestJs, typeAliases)
    val irValue =
      ValueSerde.deserializeWithType(manifestJs.fields(Archive.ManifestValueKey),
                                     irType,
                                     Archive.ManifestValueKey)
    (irType, irValue)
  }

  override lazy val (irType: Type, irValue: Value) = {
    packedTypeAndValue.getOrElse(readManifest)
  }

  /**
    * Unpacks the files in the archive relative to the given parent dir, and updates
    * the paths within `irValue` and returns a new Archive object.
    * @param mount whether to mount the filesystem
    * @param name the name of the variable associated with this archive
    * @return the updated Archive object and a Vector of localized paths
    */
  def localize(mount: Boolean = true,
               name: Option[String] = None): (LocalizedArchive, Vector[Path]) = {
    def transformer(relPath: Path): Path = {
      archive.resolve(relPath)
    }

    // make sure these lazy vals have been instantiated
    val (t, v) = (irType, irValue)
    val (localizedValue, filePaths) = Archive.transformPaths(v, t, transformer)
    if (mount && !archive.isMounted) {
      archive.mount()
    }
    val localizedArchive =
      LocalizedArchive(t, localizedValue)(Some(archive, v), Some(archive.mountPoint), name)
    (localizedArchive, filePaths.values.toVector)
  }

  def isOpen: Boolean = {
    archive.isMounted
  }

  def close(): Unit = {
    archive.unmount()
  }
}

/**
  * Represents an archive file that has been localized, i.e. the files referenced
  * in its `irValue` have absolute paths to files within a squashfs image that is
  * mounted on disk. It could be a previously packed archive or a new archive.
  * @param irType WdlType
  * @param irValue WdlValue
  * @param encoding character encoding of contained files
  * @param packedArchiveAndValue optional archive file and delocalized value of a
  *                              PackedArchive from which this LocalizedArchive
  *                              was created
  * @param parentDir optional Path to which files in the archive are relativeized
  *                  when packing.
  * @param name an optional name that will be used to prefix the randomly-generated
  *             archive name, if `originalPath` is `None`
  */
case class LocalizedArchive(
    irType: Type,
    irValue: Value,
    encoding: Charset = FileUtils.DefaultEncoding
)(packedArchiveAndValue: Option[(SquashFs, Value)] = None,
  parentDir: Option[Path] = None,
  name: Option[String] = None)
    extends Archive {
  assert(packedArchiveAndValue.isDefined || parentDir.isDefined)
  override val localized: Boolean = true

  private lazy val archive: SquashFs = {
    packedArchiveAndValue
      .map(_._1)
      .getOrElse(
          SquashFs(
              Files.createTempFile(name.getOrElse(Archive.ArchiveFilePrefix),
                                   Archive.ArchiveFileSuffix)
          )()
      )
  }

  override lazy val path: Path = archive.image

  private def createArchive(t: Type,
                            v: Value,
                            files: Map[Path, Path],
                            removeSourceFiles: Boolean): Unit = {
    val manifest = JsObject(
        TypeSerde.serializeOne(t).fields ++ Map(Archive.ManifestValueKey -> ValueSerde.serialize(v))
    )
    try {
      // write the manifest to a temp file
      val manifestDir = Files.createTempDirectory("manifest")
      val manifestPath = manifestDir.resolve(Archive.ManifestFile)
      JsUtils.jsToFile(manifest, manifestPath)
      // add the manifest to the root of the archive
      archive.append(manifestPath, removeSource = removeSourceFiles)
      // write each group of files, where files are grouped
      // by their target subdirectory within the archive
      files.groupMap(_._2.getParent)(_._1).foreach {
        case (relParent, srcFiles) =>
          archive.appendAll(srcFiles.toVector, Option(relParent), removeSource = removeSourceFiles)
      }
    } catch {
      case ex: Throwable =>
        throw new RuntimeException(s"error writing archive to ${path}", ex)
    }
  }

  def pack(removeSourceFiles: Boolean = true): PackedArchive = {
    val delocalizedValue = packedArchiveAndValue.map(_._2).getOrElse {
      def transformer(absPath: Path): Path = {
        parentDir
          .map { parent =>
            if (!absPath.startsWith(parent)) {
              throw new RuntimeException(
                  s"path ${absPath} is not located under parent dir ${parentDir}"
              )
            }
            parent.relativize(absPath)
          }
          .getOrElse(throw new RuntimeException("parentDir is required to relativize local paths"))
      }
      val (delocalizedValue, filePaths) = Archive.transformPaths(irValue, irType, transformer)
      createArchive(irType, delocalizedValue, filePaths, removeSourceFiles)
      delocalizedValue
    }
    PackedArchive(path, encoding)(packedTypeAndValue = Some(irType, delocalizedValue))
  }

  def isOpen: Boolean = {
    archive.isMounted
  }

  def close(): Unit = {
    if (archive.isMounted) {
      archive.unmount()
    }
  }

  def toJson: JsValue = {
    JsObject(
        "type" -> TypeSerde.serializeOne(irType),
        "value" -> ValueSerde.serialize(irValue),
        "packedValue" -> packedArchiveAndValue
          .map {
            case (_, v) => ValueSerde.serialize(v)
          }
          .getOrElse(JsNull),
        "fs" -> archive.toJson,
        "encoding" -> JsString(encoding.name()),
        "parentDir" -> parentDir.map(p => JsString(p.toString)).getOrElse(JsNull),
        "name" -> name.map(JsString(_)).getOrElse(JsNull)
    )
  }
}

object LocalizedArchive {
  def fromJson(jsValue: JsValue, schemas: Map[String, TSchema]): LocalizedArchive = {
    jsValue.asJsObject
      .getFields("type", "value", "packedValue", "fs", "encoding", "parentDir", "name") match {
      case Seq(typeJs, valueJs, packedValueJs, fsJs, encodingJs, parentDirJs, nameJs) =>
        val (irType, _) = TypeSerde.deserializeOne(typeJs, schemas)
        val irValue = ValueSerde.deserializeWithType(valueJs, irType)
        val archive = SquashFs.fromJson(fsJs)
        val packedArchiveAndValue = packedValueJs match {
          case JsNull => None
          case v =>
            val value = ValueSerde.deserializeWithType(v, irType)
            Some((archive, value))
        }
        val encoding = Charset.forName(JsUtils.getString(encodingJs))
        val parentDir = parentDirJs match {
          case JsString(path) => Some(Paths.get(path))
          case JsNull         => None
          case _              => throw new Exception(s"invalid parentDir ${parentDirJs}")
        }
        val name = nameJs match {
          case JsString(name) => Some(name)
          case JsNull         => None
          case _              => throw new Exception(s"invalid name ${nameJs}")
        }
        LocalizedArchive(irType, irValue, encoding)(packedArchiveAndValue, parentDir, name)
      case _ =>
        throw new Exception(s"invalid serialized LocalizedArchive value ${jsValue}")
    }
  }
}
