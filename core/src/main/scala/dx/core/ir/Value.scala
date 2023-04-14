package dx.core.ir

import dx.core.ir.Type._
import dx.util.CollectionUtils.IterableOnceExtensions
import dx.util.protocols.DxFolderSource

import scala.collection.immutable.SeqMap

/**
  * A language-independent representation of values used as input to/output from
  * DNAnexus applications and workflows.
  */
sealed trait Value

object Value {
  // Primitive values.
  sealed trait PrimitiveValue extends Value
  case class VInt(value: Long) extends PrimitiveValue
  case class VFloat(value: Double) extends PrimitiveValue
  case class VString(value: String) extends PrimitiveValue
  case class VBoolean(value: Boolean) extends PrimitiveValue

  sealed trait PathValue extends PrimitiveValue

  /**
    * Represents a local or remote file. If `contents` is set,
    * this is a virtual file can be localized by writing its
    * contents to a local path.
    * @param uri the local file path or remote URI
    * @param basename the name to use when localizing the file
    * @param contents the file contents
    * @param checksum the file checksum
    * @param size the file size
    * @param secondaryFiles additional files that must accompany this file
    * @param format File format according to the language specification
    *
    */
  case class VFile(uri: String,
                   basename: Option[String] = None,
                   contents: Option[String] = None,
                   checksum: Option[String] = None,
                   size: Option[Long] = None,
                   secondaryFiles: Vector[PathValue] = Vector.empty,
                   format: Option[String] = None,
                   metadata: Option[String] = None)
      extends PathValue

  sealed trait DirectoryValue extends PathValue

  /**
    * A local directory or remote DNAnexus folder.
    * @param uri local path or dx://project-xxx:/path/to/folder/ URI
    * @param basename the name to use when localizing the directory
    * @param listing an Optional Vector of files/subfolders in this folder.
    *              If specified, only these files/subfolders will be
    *              localized.
    */
  case class VFolder(uri: String,
                     basename: Option[String] = None,
                     listing: Option[Vector[PathValue]] = None)
      extends DirectoryValue

  /**
    * A synthetic directory consisting of specific platform files/folders.
    * @param basename the name to use when localizing the directory
    * @param items a Vector of files/subdirectories in this directory
    */
  case class VListing(basename: String, items: Vector[PathValue] = Vector.empty)
      extends DirectoryValue

  /**
    * Represents the empty value for an optional field.
    */
  case object VNull extends Value

  /**
    * Represents the empty value for an optional field which has been forced in frag input evaluation.
    */
  case object VForcedNull extends Value

  /**
    * An array of values
    */
  case class VArray(items: Vector[Value]) extends Value

  object VArray {
    val empty: VArray = VArray(Vector())

    def apply(items: Value*): VArray = {
      VArray(items.toVector)
    }
  }

  /**
    * A JSON object. Fields are stored in a SeqMap to preserve their order.
    */
  case class VHash(fields: SeqMap[String, Value]) extends Value

  object VHash {
    def apply(fields: (String, Value)*): VHash = {
      new VHash(fields.to(SeqMap))
    }
  }

  trait WalkHandler[T] {
    def apply(value: Value, t: Option[Type], optional: Boolean, ctx: T): Option[T]
  }

  def walk[T](value: Value, t: Option[Type], initialContext: T, handler: WalkHandler[T]): T = {
    def walkPaths(paths: Vector[PathValue], innerContext: T): T = {
      paths
        .foldLeft(innerContext) {
          case (ctx, f: VFile)          => inner(f, Some(TFile), ctx)
          case (ctx, d: DirectoryValue) => inner(d, Some(TDirectory), ctx)
        }
    }
    def inner(innerValue: Value, innerType: Option[Type], innerContext: T): T = {
      val (nonOptType, optional) = innerType match {
        case Some(t) if isOptional(t) => (Some(unwrapOptional(innerType.get)), true)
        case Some(_)                  => (innerType, false)
        case None                     => (innerType, true)
      }
      handler(innerValue, nonOptType, optional, innerContext).getOrElse {
        (nonOptType, innerValue) match {
          case (_, VNull) if optional       => innerContext
          case (_, VForcedNull) if optional => innerContext
          case (_, VNull) =>
            throw new Exception(s"null value for non-optional type ${innerType.get}")
          case (_, VForcedNull) =>
            throw new Exception(s"null value for non-optional type ${innerType.get}")
          case (Some(TFile) | None, f: VFile) if f.secondaryFiles.nonEmpty =>
            walkPaths(f.secondaryFiles, innerContext)
          case (Some(TDirectory) | None, f: VFolder) if f.listing.exists(_.nonEmpty) =>
            walkPaths(f.listing.get, innerContext)
          case (Some(TDirectory) | None, l: VListing) if l.items.nonEmpty =>
            walkPaths(l.items, innerContext)
          case (Some(TArray(_, true)), VArray(Vector())) =>
            throw new Exception("empty array for non-empty array type")
          case (Some(TArray(itemType, _)), VArray(items)) =>
            items.foldLeft(innerContext) {
              case (ctx, item) => inner(item, Some(itemType), ctx)
            }
          case (_, VArray(items)) =>
            items.foldLeft(innerContext) {
              case (ctx, item) => inner(item, None, ctx)
            }
          case (Some(TSchema(name, fieldTypes)), VHash(fields)) =>
            fields.foldLeft(innerContext) {
              case (_, (k, _)) if !fieldTypes.contains(k) =>
                throw new Exception(s"invalid member ${k} of schema ${name}")
              case (ctx, (k, v)) => inner(v, Some(fieldTypes(k)), ctx)
            }
          case (_, VHash(fields)) =>
            fields.foldLeft(innerContext) {
              case (ctx, (_, v)) => inner(v, None, ctx)
            }
          case (Some(TEnum(symbols)), s: VString) if symbols.contains(s.value) => innerContext
          case (Some(TEnum(symbols)), other) =>
            throw new Exception(
                s"${other} is not one of the allowed symbols ${symbols.mkString(",")}"
            )
          case (Some(TMulti.Any), _) => innerContext
          case (Some(TMulti(bounds)), _) =>
            bounds
              .collectFirstDefined { boundType =>
                try {
                  Some(inner(innerValue, Some(boundType), innerContext))
                } catch {
                  case _: Throwable => None
                }
              }
              .getOrElse(
                  throw new Exception(s"could not handle ${innerValue} as any of ${bounds}")
              )
          case _ => innerContext
        }
      }
    }
    inner(value, t, initialContext)
  }

  trait TransformHandler {
    def apply(value: Value, t: Option[Type], optional: Boolean): Option[Value]
  }

  /**
    * Transforms a Value to another Value, applying the `handler` function
    * at each level of nesting. The default rules handle recursively
    * descending into parameterized types (VArray, VHash) but otherwise
    * leave the original value unchanged.
    * @param value the Value to transform
    * @param t an optional Type to which the value should be transformed
    * @param handler a function that may transform a Value to another Value.
    *                If the function returns Some(newValue), then newValue is
    *                the result of the transformation, otherwise the default
    *                transformation rules are applied.
    * @return the transformed Value
    */
  def transform(value: Value, t: Option[Type], handler: TransformHandler): Value = {
    def transformPaths(paths: Vector[PathValue]): Vector[PathValue] = {
      paths
        .map {
          case f: VFile          => inner(f, Some(TFile))
          case d: DirectoryValue => inner(d, Some(TDirectory))
        }
        .map {
          case p: PathValue => p
          case other        => throw new Exception(s"not a PathValue ${other}")
        }
    }
    def inner(innerValue: Value, innerType: Option[Type] = None): Value = {
      val (nonOptType, optional) = innerType match {
        case Some(t) if isOptional(t) => (Some(unwrapOptional(innerType.get)), true)
        case Some(_)                  => (innerType, false)
        case None                     => (innerType, true)
      }
      handler(innerValue, nonOptType, optional).getOrElse {
        (nonOptType, innerValue) match {
          case (_, VNull) if optional       => VNull
          case (_, VForcedNull) if optional => VForcedNull
          case (_, VNull) =>
            throw new Exception(s"null value for non-optional type ${innerType.get}")
          case (_, VForcedNull) =>
            throw new Exception(s"null value for non-optional type ${innerType.get}")
          case (Some(TFile) | None, f: VFile) if f.secondaryFiles.nonEmpty =>
            f.copy(secondaryFiles = transformPaths(f.secondaryFiles))
          case (Some(TDirectory) | None, f: VFolder) if f.listing.exists(_.nonEmpty) =>
            f.copy(listing = Some(transformPaths(f.listing.get)))
          case (Some(TDirectory) | None, l: VListing) if l.items.nonEmpty =>
            l.copy(items = transformPaths(l.items))
          case (Some(TArray(_, true)), VArray(Vector())) =>
            throw new Exception("empty array for non-empty array type")
          case (Some(TArray(itemType, _)), VArray(items)) =>
            VArray(items.map(inner(_, Some(itemType))))
          case (_, VArray(items)) => VArray(items.map(inner(_)))
          case (Some(TSchema(name, fieldTypes)), VHash(fields)) =>
            VHash(fields.map {
              case (k, _) if !fieldTypes.contains(k) =>
                throw new Exception(s"invalid member ${k} of schema ${name}")
              case (k, v) => k -> inner(v, Some(fieldTypes(k)))
            })
          case (_, VHash(fields)) =>
            VHash(fields.map { case (k, v) => k -> inner(v) })
          case (Some(TEnum(symbols)), s: VString) if symbols.contains(s.value) => s
          case (Some(TEnum(symbols)), other) =>
            throw new Exception(
                s"${other} is not one of the allowed symbols ${symbols.mkString(",")}"
            )
          case (Some(TMulti.Any), _) => innerValue
          case (Some(TMulti(bounds)), _) =>
            bounds
              .collectFirstDefined { boundType =>
                try {
                  Some(inner(innerValue, Some(boundType)))
                } catch {
                  case _: Throwable => None
                }
              }
              .getOrElse(
                  throw new Exception(s"could not transform ${innerValue} as any of ${bounds}")
              )
          case _ => innerValue
        }
      }
    }
    inner(value, t)
  }

  def coerceTo(value: Value, targetType: Type): Value = {
    (targetType, value) match {
      case (TOptional(_), VNull)       => VNull
      case (TOptional(_), VForcedNull) => VForcedNull
      case (TOptional(t), _)           => coerceTo(value, t)
      // check whether the value is already of the correct type
      case (TBoolean, b: VBoolean)    => b
      case (TInt, i: VInt)            => i
      case (TFloat, f: VFloat)        => f
      case (TString, s: VString)      => s
      case (TFile, f: VFile)          => f
      case (TDirectory, p: PathValue) => p
      case (THash, h: VHash)          => h
      // compound types
      case (TArray(_, nonEmpty), VArray(items)) if nonEmpty && items.isEmpty =>
        throw new Exception("cannot coerce empty array to non-empty array")
      case (TArray(t, _), VArray(items)) =>
        VArray(items.map(coerceTo(_, t)))
      case (TSchema(schemaName, fields), VHash(members)) =>
        val invalid = members.keySet.diff(fields.keySet)
        if (invalid.nonEmpty) {
          throw new Exception(
              s"cannot coerce hash with key(s) ${invalid.mkString(",")} to type ${targetType}"
          )
        }
        VHash(fields.collect {
          case (fieldName, t) if members.contains(fieldName) =>
            fieldName -> coerceTo(members(fieldName), t)
          case (fieldName, t) if !Type.isOptional(t) =>
            throw new Exception(s"missing required ${schemaName} field ${fieldName}")
        })
      case (TEnum(symbols), s: VString) if symbols.contains(s.value) => s
      case (TEnum(symbols), _) =>
        throw new Exception(s"${value} is not one of allowed symbols ${symbols.mkString(",")}")
      case (TMulti.Any, _) => value
      case (TMulti(bounds), _) =>
        bounds.iterator
          .collectFirstDefined { t =>
            try {
              Some(coerceTo(value, t))
            } catch {
              case _: Throwable => None
            }
          }
          .getOrElse(
              throw new Exception(s"cannot coerce ${value} to any of ${bounds.mkString(",")}")
          )
      // coercions
      case (TFile, VString(s))                                         => VFile(s)
      case (TString, f: VFile)                                         => VString(f.uri)
      case (TDirectory, VString(s)) if DxFolderSource.isDxFolderUri(s) => VFolder(s)
      case (TString, f: VFolder)                                       => VString(f.uri)
      case (TFloat, VInt(i))                                           => VFloat(i.toFloat)
      case _ =>
        throw new Exception(s"cannot coerce ${value} to ${targetType}")
    }
  }
}
