package dx.core.ir

import java.nio.charset.Charset
import java.nio.file.{Files, Path, Paths}

import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.{FileUtils, JsUtils, Logger, SysUtils}
import spray.json._

case class SquashFs(archiveFile: Path)(
    val mountPoint: Path = Files.createTempDirectory(archiveFile.getFileName.toString),
    logger: Logger = Logger.get
) {
  private var mounted: Boolean = false
  private lazy val isLinux = System.getProperty("os.name").toLowerCase.startsWith("linux")

  def isMounted: Boolean = mounted

  /**
    * Mounts the archive if it has not already been mounted.
    */
  def mount(): Unit = {
    if (!mounted) {
      try {
        if (isLinux) {
          SysUtils.execCommand(s"mount -t squashfs ${archiveFile.toString} ${mountPoint.toString}")
        } else {
          SysUtils.execCommand(s"squashfuse ${archiveFile.toString} ${mountPoint.toString}")
        }
        mounted = true
      } catch {
        case ex: Throwable =>
          throw new RuntimeException(s"unable to mount archive ${archiveFile} at ${mountPoint}", ex)
      }
    }
  }

  /**
    * Unmounts the archive if it has previously been mounted.
    */
  def unmount(): Unit = {
    if (isMounted) {
      try {
        SysUtils.execCommand(s"umount ${mountPoint.toString}")
        mounted = false
      } catch {
        case ex: Throwable =>
          logger.error(s"error unmounting archive at ${mountPoint}", Some(ex))
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
    * @param path path to the file to add to the archive
    * @param parentDir parent directory to which to relativize `path` - if None, the
    *                  file is added to the root of the fs
    * @param removeOriginal remove the original file after it has been added to the
    *                       archive
    */
  def append(path: Path, parentDir: Option[Path] = None, removeOriginal: Boolean = false): Unit = {
    try {
      if (parentDir.isDefined && parentDir != path.getParent) {
        // mksquashfs can't easily add a specific file in a subdir, so we need to
        // create a temp directory with the desired sub-directory structure, then
        // create a hard link to the file within the leaf dir, then append the file,
        // and finally delete the link
        val lnDir = appendDir.resolve(parentDir.get.relativize(path.getParent))
        val lnFile = lnDir.resolve(path.getFileName)
        SysUtils.execCommand(
            s"""mkdir -p ${lnDir.toString}
               |&& ln ${path} ${lnFile.toString}
               |&& mksquashfs ${appendDir.toString} ${archiveFile.toString} -quiet
               |&& rm -Rf ${appendDir}/*""".stripMargin.replaceAll("\n", " ")
        )
      } else {
        SysUtils.execCommand(s"mksquashfs ${path} ${archiveFile} -quiet")
      }
    } catch {
      case ex: Throwable =>
        throw new RuntimeException(s"error adding ${path} to archive ${archiveFile}", ex)
    }
    if (removeOriginal) {
      try {
        Files.delete(path)
      } catch {
        case ex: Throwable =>
          logger.error(s"error deleting ${path} after adding it to archive ${archiveFile}",
                       Some(ex))
      }
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
  val ManifestFile: String = "manifest.json"
  val ManifestValueKey: String = "value"

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
      (VFile(newPath.toString), Map(oldPath -> newPath))
    }
    def transformNoType(innerValue: Value): (Value, Map[Path, Path]) = {
      innerValue match {
        case VFile(s) => transformFile(s)
        case VArray(value) =>
          val (items, paths) = value.map(transformNoType).unzip
          (VArray(items), paths.flatten.toMap)
        case VHash(members) =>
          val (newMembers, paths) = members.map {
            case (k, v) =>
              val (value, paths) = transformNoType(v)
              (k -> value, paths)
          }.unzip
          (VHash(newMembers.toMap), paths.flatten.toMap)
        case _: VArchive =>
          throw new RuntimeException("nested archive values are not allowed")
        case _ =>
          (innerValue, Map.empty)
      }
    }
    def transformWithType(innerValue: Value, innerType: Type): (Value, Map[Path, Path]) = {
      (innerType, innerValue) match {
        case (TFile, VFile(s))   => transformFile(s)
        case (TFile, VString(s)) => transformFile(s)
        case (TOptional(t), value) if value != VNull =>
          val (v, paths) = transformWithType(value, t)
          (v, paths)
        case (TArray(itemType, _), VArray(value)) =>
          val (items, paths) = value.map(transformWithType(_, itemType)).unzip
          (VArray(items), paths.flatten.toMap)
        case (TSchema(name, memberTypes), VHash(members)) =>
          val (newMembers, paths) = members.map {
            case (k, v) =>
              val valueType = memberTypes.getOrElse(
                  k,
                  throw new RuntimeException(s"${k} is not a member of schema ${name}")
              )
              val (value, paths) = transformWithType(v, valueType)
              (k -> value, paths)
          }.unzip
          (VHash(newMembers.toMap), paths.flatten.toMap)
        case (THash, _: VHash) =>
          transformNoType(innerValue)
        case (_, _: VArchive) =>
          throw new RuntimeException("nested archive values are not allowed")
        case _ =>
          (innerValue, Map.empty)
      }
    }
    transformWithType(irValue, irType)
  }
}

case class PackedArchive(path: Path, encoding: Charset = FileUtils.DefaultEncoding)(
    typeAliases: Option[Map[String, TSchema]] = None,
    packedTypeAndValue: Option[(Type, Value)] = None
) extends Archive {
  assert(typeAliases.isDefined || packedTypeAndValue.isDefined)
  if (!Files.exists(path)) {
    throw new Exception(s"${path} does not exist")
  } else if (Files.isDirectory(path)) {
    throw new Exception(s"${path} is not a file")
  }

  override val localized: Boolean = false
  private val archive: SquashFs = SquashFs(path)()
  private lazy val manifestPath: Path = archive.resolve(Archive.ManifestFile)

  private def readManifest: (Type, Value) = {
    if (!archive.isMounted) {
      archive.mount()
    }
    val manifestJs =
      try {
        JsUtils.jsFromFile(manifestPath).asJsObject
      } catch {
        case ex: Throwable =>
          throw new RuntimeException(s"unable to read archive manifest ${manifestPath}", ex)
      }
    val (irType, _) = TypeSerde.deserializeOne(manifestJs, typeAliases.get)
    val irValue =
      ValueSerde.deserializeWithType(manifestJs.fields(Archive.ManifestValueKey), irType)
    (irType, irValue)
  }

  override lazy val (irType: Type, irValue: Value) = {
    packedTypeAndValue.getOrElse(readManifest)
  }

  /**
    * Unpacks the files in the archive relative to the given parent dir, and updates
    * the paths within `irValue` and returns a new Archive object.
    * @param name the name of the variable associated with this archive
    * @return the updated Archive object and a Vector of localized paths
    */
  def localize(name: Option[String] = None): (LocalizedArchive, Vector[Path]) = {
    def transformer(relPath: Path): Path = {
      archive.resolve(relPath)
    }

    // make sure these lazy vals have been instantiated
    val (t, v) = (irType, irValue)
    val (localizedValue, filePaths) = Archive.transformPaths(v, t, transformer)
    val unpackedArchive =
      LocalizedArchive(t, localizedValue)(Some(archive, v), Some(archive.mountPoint), name)
    (unpackedArchive, filePaths.values.toVector)
  }

  def close(): Unit = {
    archive.unmount()
  }
}

/**
  * Represents an archive file that has been localized, i.e. the files referenced
  * in its `irValue` have absolute paths to files within a squasfs image that is
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
  override val localized: Boolean = true

  private lazy val archive: SquashFs = {
    packedArchiveAndValue
      .map(_._1)
      .getOrElse(SquashFs(Files.createTempFile(name.getOrElse("archive"), ".img"))())
  }

  override lazy val path: Path = archive.archiveFile

  private def createArchive(t: Type, v: Value, filePaths: Map[Path, Path]): Unit = {
    val manifest = JsObject(
        TypeSerde.serializeOne(t).fields ++ Map(Archive.ManifestValueKey -> ValueSerde.serialize(v))
    )
    try {
      // write the manifest to a temp file
      val manifestDir = Files.createTempDirectory("manifest")
      val manifestPath = manifestDir.resolve(Archive.ManifestFile)
      JsUtils.jsToFile(manifest, manifestPath)
      // add the manifest to the archive
      archive.append(manifestPath, Some(manifestDir), removeOriginal = true)
      // write each file
      filePaths.foreach {
        case (absPath, relPath) =>
          archive.append(absPath, Some(relPath), removeOriginal = true)
      }
    } catch {
      case ex: Throwable =>
        throw new RuntimeException(s"error writing archive to ${path}", ex)
    }
  }

  lazy val pack: PackedArchive = {
    val delocalizedValue = packedArchiveAndValue.map(_._2).getOrElse {
      def transformer(absPath: Path): Path = {
        parentDir.get.relativize(absPath)
      }
      val (delocalizedValue, filePaths) = Archive.transformPaths(irValue, irType, transformer)
      createArchive(irType, delocalizedValue, filePaths)
      delocalizedValue
    }
    PackedArchive(path, encoding)(packedTypeAndValue = Some(irType, delocalizedValue))
  }
}
