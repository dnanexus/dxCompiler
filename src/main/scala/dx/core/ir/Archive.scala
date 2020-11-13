package dx.core.ir

import java.io.{BufferedInputStream, FileInputStream, FileOutputStream}
import java.nio.charset.Charset
import java.nio.file.{Files, Path, Paths}

import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.FileUtils
import org.apache.commons.compress.archivers.tar.{
  TarArchiveEntry,
  TarArchiveInputStream,
  TarArchiveOutputStream
}
import org.apache.commons.compress.utils.IOUtils
import spray.json._

/**
  * An Archive is a TAR file that contains 1) a JSON file (called manifest.json)
  * with the serialized representation of a complex value, along with its serialized
  * type information, and 2) the files referenced by all of the VFile values nested
  * within the complex value.
  * TODO: add option to support compressing the TAR
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
  private var isOpen = false

  private lazy val tarStream: TarArchiveInputStream =
    try {
      val inputStream = new FileInputStream(path.toFile)
      val tarStream = new TarArchiveInputStream(if (inputStream.markSupported()) {
        inputStream
      } else {
        new BufferedInputStream(inputStream)
      })
      isOpen = true
      tarStream
    } catch {
      case ex: Throwable =>
        throw new Exception(s"invalid WDL value archive ${path}", ex)
    }

  private object iterator extends Iterator[TarArchiveEntry] {
    private var currentEntry: TarArchiveEntry = _

    override def hasNext: Boolean = {
      currentEntry = tarStream.getNextTarEntry
      currentEntry != null
    }

    override def next(): TarArchiveEntry = currentEntry
  }

  private def readManifest: (Type, Value) = {
    if (isOpen) {
      throw new RuntimeException("manifest has already been read")
    }
    if (!iterator.hasNext) {
      throw new RuntimeException(s"invalid archive file: missing ${Archive.ManifestFile}")
    }
    val manifestEntry = iterator.next()
    if (manifestEntry.isFile && manifestEntry.getName == Archive.ManifestFile) {
      val contents = new String(IOUtils.toByteArray(tarStream), encoding)
      val manifestJs = contents.parseJson.asJsObject
      val (irType, _) = TypeSerde.deserializeOne(manifestJs, typeAliases.get)
      val irValue =
        ValueSerde.deserializeWithType(manifestJs.fields(Archive.ManifestValueKey), irType)
      (irType, irValue)
    } else {
      throw new RuntimeException(
          s"invalid archive file: expected first entry to be ${Archive.ManifestFile}, not ${manifestEntry}"
      )
    }
  }

  override lazy val (irType: Type, irValue: Value) = {
    packedTypeAndValue.getOrElse(readManifest)
  }

  /**
    * Unpacks the files in the archive relative to the given parent dir, and updates
    * the paths within `irValue` and returns a new Archive object.
    * @param parentDir the directory in which to localize files
    * @return the updated Archive object and a Vector of localized paths
    */
  def localize(parentDir: Path, name: Option[String] = None): (LocalizedArchive, Vector[Path]) = {
    def transformer(relPath: Path): Path = {
      parentDir.resolve(relPath)
    }

    // make sure these lazy vals have been instantiated
    val (t, v) = (irType, irValue)
    val (localizedValue, filePaths) = Archive.transformPaths(v, t, transformer)
    val unpackedArchive = LocalizedArchive(t, localizedValue)(Some(path, v), Some(parentDir), name)
    (unpackedArchive, filePaths.values.toVector)
  }

  def close(): Unit = {
    if (isOpen) {
      tarStream.close()
    }
  }
}

/**
  * Represents an archive file that has been localized, i.e. the files referenced
  * in its `irValue` are localized on disk. It could be a previously packed
  * archive or a new archive that has not yet been packed.
  *
  * Either `packedPathAndValue` or `parentDir` must be specified.
  * @param irType WdlType
  * @param irValue WdlValue
  * @param encoding character encoding of contained files
  * @param packedPathAndValue optional path and delocalized value of a PackedArchive
  *                           from which this LocalizedArchive was created
  * @param parentDir optional Path to which files in the archive are relativeized
  *                  when packing.
  * @param name an optional name that will be used to prefix the randomly-generated
  *             archive name, if `originalPath` is `None`
  */
case class LocalizedArchive(
    irType: Type,
    irValue: Value,
    encoding: Charset = FileUtils.DefaultEncoding
)(packedPathAndValue: Option[(Path, Value)] = None,
  parentDir: Option[Path] = None,
  name: Option[String] = None)
    extends Archive {
  assert(packedPathAndValue.isDefined || parentDir.isDefined)
  override val localized: Boolean = true

  override lazy val path: Path = {
    packedPathAndValue
      .map(_._1)
      .getOrElse(Files.createTempFile(name.getOrElse("archive"), ".tar"))
  }

  private def createArchive(path: Path, t: Type, v: Value, filePaths: Map[Path, Path]): Unit = {
    val manifest = JsObject(
        TypeSerde.serializeOne(t).fields ++ Map(Archive.ManifestValueKey -> ValueSerde.serialize(v))
    )
    val manifestBytes = manifest.prettyPrint.getBytes(encoding)
    val tarStream = new TarArchiveOutputStream(new FileOutputStream(path.toFile))
    try {
      // write the manifest
      val manifestEntry = new TarArchiveEntry(Archive.ManifestFile)
      manifestEntry.setSize(manifestBytes.size)
      tarStream.putArchiveEntry(manifestEntry)
      tarStream.write(manifestBytes)
      tarStream.closeArchiveEntry()
      // write each file
      filePaths.foreach {
        case (absPath, relPath) =>
          val dirEntry = new TarArchiveEntry(absPath.toFile, relPath.toString)
          tarStream.putArchiveEntry(dirEntry)
          Files.copy(absPath, tarStream)
          tarStream.closeArchiveEntry()
      }
    } finally {
      tarStream.close()
    }
  }

  lazy val pack: PackedArchive = {
    val delocalizedValue = packedPathAndValue.map(_._2).getOrElse {
      def transformer(absPath: Path): Path = {
        parentDir.get.relativize(absPath)
      }
      val (delocalizedValue, filePaths) = Archive.transformPaths(irValue, irType, transformer)
      createArchive(path, irType, delocalizedValue, filePaths)
      delocalizedValue
    }
    PackedArchive(path, encoding)(packedTypeAndValue = Some(irType, delocalizedValue))
  }
}
