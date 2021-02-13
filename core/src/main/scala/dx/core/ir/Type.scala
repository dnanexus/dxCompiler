package dx.core.ir

import scala.annotation.tailrec
import scala.collection.immutable.SeqMap

sealed trait Type

object Type {
  // Primitive types that are supported natively.
  sealed trait PrimitiveType extends Type
  case object TBoolean extends PrimitiveType
  case object TInt extends PrimitiveType
  case object TFloat extends PrimitiveType
  case object TString extends PrimitiveType
  case object TFile extends PrimitiveType

  /**
    * A directory maps to a DNAnexus folder, or to some other representation of
    * a hierarchy of files, e.g. a tar or zip archive.
    */
  case object TDirectory extends PrimitiveType

  /**
    * Wrapper that indicates a type is optional.
    * @param t wrapped type
    */
  case class TOptional(t: Type) extends Type

  sealed trait TCollection extends Type

  /**
    * Array of primitive or file values - all items in an array must be of the
    * same type. Some languages (e.g. WDL) have a quantifier to specify that
    * the array must be non-empty - this does not change how the array is
    * represented, but an error may be thrown if the array value is empty.
    * @param t inner type
    * @param nonEmpty whether the array must not be empty
    */
  case class TArray(t: Type, nonEmpty: Boolean = false) extends TCollection

  /**
    * A JSON object.
    */
  case object THash extends TCollection

  /**
    * Represents a user-defined type. Values of this type are represented
    * as VHash. Fields are stored in a SeqMap to preserve their order.
    * @param name type name
    */
  case class TSchema(name: String, fields: SeqMap[String, Type]) extends TCollection

  /**
    * Represents a String that is restricted to the specified symbols.
    * Symbols are stored in a Vector to preserve their order.
    * @param symbols set of allowed symbols in enum
    */
  case class TEnum(symbols: Vector[String]) extends PrimitiveType

  /**
    * Represents a collection of equally valid types. If bounds
    * is empty, then any type is allowed.
    */
  case class TMulti(bounds: Vector[Type]) extends Type

  object TMulti {
    val Any: TMulti = TMulti(Vector())
  }

  @tailrec
  def isPrimitive(t: Type): Boolean = {
    t match {
      case _: PrimitiveType => true
      case TOptional(inner) => isPrimitive(inner)
      case _                => false
    }
  }

  def isNativePrimitive(t: Type): Boolean = {
    t match {
      case TBoolean => true
      case TInt     => true
      case TFloat   => true
      case TString  => true
      case TFile    => true
      // TODO: TDirectory
      case _ => false
    }
  }

  /**
    * Is this an IR type that maps to a native DX type?
    * DNAnexus only supports primitives, optionals of primitives, and arrays of primitives
    * (with no nested optional types).
    * @param t IR type
    * @return
    */
  def isNative(t: Type): Boolean = {
    t match {
      case TOptional(TArray(inner, _)) =>
        // TODO: should an optional non-empty array be considered native? Currently it is
        //  allowed, and must always agree with the type conversion logic in TypeSerde.toNative.
        isNativePrimitive(inner)
      case TOptional(inner) => isNativePrimitive(inner)
      case TArray(inner, _) => isNativePrimitive(inner)
      case _                => isNativePrimitive(t)
    }
  }

  def isOptional(t: Type): Boolean = {
    t match {
      case _: TOptional  => true
      case TMulti(types) =>
        // a multi-type is considered optional if any of its
        // alternative types is optional, regardless of whether
        // it is wrapped in TOptional
        types.exists(isOptional)
      case _ => false
    }
  }

  /**
    * Makes a type optional.
    * @param t the type
    * @param force if true, then `t` will be made optional even if it is already optional.
    * @return
    */
  def ensureOptional(t: Type, force: Boolean = false): TOptional = {
    t match {
      case t if force   => TOptional(t)
      case t: TOptional => t
      case _            => TOptional(t)
    }
  }

  def unwrapOptional(t: Type, mustBeOptional: Boolean = false): Type = {
    t match {
      case TOptional(wrapped) => wrapped
      case _ if mustBeOptional && !isOptional(t) =>
        throw new Exception(s"Type ${t} is not TOptional")
      case _ => t
    }
  }

  def isNestedOptional(t: Type): Boolean = {
    t match {
      case TOptional(TOptional(_)) => true
      case _                       => false
    }
  }

  def collectSchemas(types: Vector[Type]): Map[String, TSchema] = {
    def inner(t: Type, schemas: Map[String, TSchema]): Map[String, TSchema] = {
      t match {
        case schema: TSchema if !schemas.contains(schema.name) =>
          schema.fields.values.foldLeft(schemas + (schema.name -> schema)) {
            case (accu, memberType) => inner(memberType, accu)
          }
        case TOptional(c: TCollection) => inner(c, schemas)
        case TArray(c: TCollection, _) => inner(c, schemas)
        case _                         => schemas
      }
    }
    types.foldLeft(Map.empty[String, TSchema]) {
      case (accu, wdlType) => inner(wdlType, accu)
    }
  }
}
