package dx.core.ir

import dx.core.ir.Type._
import dx.util.CollectionUtils.IterableOnceExtensions
import scala.collection.immutable.{SeqMap, TreeSeqMap}

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
  case class VFile(value: String) extends PrimitiveValue

  /**
    * Represents a DNAnexus folder or a file-based representation of a directory
    * (e.g. a zip or tar archive file).
    * @param value directory
    */
  case class VDirectory(value: String) extends PrimitiveValue

  /**
    * Represents the empty value for an optional field.
    */
  case object VNull extends Value

  /**
    * An array of values
    */
  case class VArray(items: Vector[Value]) extends Value

  object VArray {
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
      new VHash(fields.to(TreeSeqMap))
    }
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
  def transform(value: Value,
                t: Option[Type],
                handler: (Value, Option[Type], Boolean) => Option[Value]): Value = {
    def inner(innerValue: Value, innerType: Option[Type] = None): Value = {
      val (nonOptType, optional) = if (innerType.exists(isOptional)) {
        (Some(unwrapOptional(innerType.get)), true)
      } else {
        (innerType, true)
      }
      handler(innerValue, nonOptType, optional).getOrElse {
        (nonOptType, innerValue) match {
          case (_, VNull) if optional => VNull
          case (_, VNull) =>
            throw new Exception(s"null value for non-optional type ${innerType.get}")
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
          case _ => innerValue
        }
      }
    }
    inner(value, t)
  }

  def coerceTo(value: Value, targetType: Type): Value = {
    (targetType, value) match {
      case (TOptional(_), VNull) => VNull
      case (TOptional(t), _)     => coerceTo(value, t)
      // check whether the value is already of the correct type
      case (TBoolean, b: VBoolean)     => b
      case (TInt, i: VInt)             => i
      case (TFloat, f: VFloat)         => f
      case (TString, s: VString)       => s
      case (TFile, f: VFile)           => f
      case (TDirectory, d: VDirectory) => d
      case (THash, h: VHash)           => h
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
      case (TFile, VString(s))      => VFile(s)
      case (TString, VFile(f))      => VString(f)
      case (TDirectory, VString(s)) => VDirectory(s)
      case (TString, VDirectory(d)) => VString(d)
      case (TFloat, VInt(i))        => VFloat(i.toFloat)
      case _ =>
        throw new Exception(s"cannot coerce ${value} to ${targetType}")
    }
  }
}
