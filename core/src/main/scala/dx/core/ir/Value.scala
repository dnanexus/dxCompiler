package dx.core.ir

import dx.core.ir.Type.{TArray, TEnum, TSchema}

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
      val (nonOptType, optional) = if (innerType.exists(Type.isOptional)) {
        (Some(Type.unwrapOptional(innerType.get)), true)
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
}
