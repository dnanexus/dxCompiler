package dx.core.ir

import dx.core.ir.Type.{TArray, TEnum, TMulti, TSchema}
import dx.core.ir.ValueSerde.ValueSerdeException

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

  /**
    * A JSON object. Fields are stored in a SeqMap to preserve their order.
    */
  case class VHash(fields: SeqMap[String, Value]) extends Value

  /**
    * Transform a Value to another Value.
    * @param value the Value to transform
    * @param t an optional Type to which the value should be transformed
    * @param handler an optional function to handle special cases. If the
    *                function returns Some(newValue), then newValue is
    *                the result of the transformation, otherwise the default
    *                transformation rules are used.
    * @return the transformed Value
    */
  def transform(value: Value,
                t: Option[Type],
                handler: (Value, Boolean) => Option[Value]): Value = {
    def inner(innerValue: Value, innerType: Option[Type] = None): Value = {
      val (nonOptType, optional) = if (innerType.exists(Type.isOptional)) {
        (Some(Type.unwrapOptional(innerType.get)), true)
      } else {
        (innerType, false)
      }
      handler(innerValue, optional).getOrElse {
        (nonOptType, innerValue) match {
          case (Some(TMulti(types)), value) if types.isEmpty => inner(value)
          case (Some(TMulti(types)), value: Value) =>
            types.foreach { t =>
              try {
                return inner(value, Some(t))
              } catch {
                case _: ValueSerdeException => ()
              }
            }
            throw ValueSerdeException(s"value ${value} does not match any of ${types}")
          case (Some(TArray(_, true)), VArray(Vector())) =>
            throw ValueSerdeException("empty array for non-empty array type")
          case (Some(TArray(itemType, _)), VArray(items)) =>
            VArray(items.map(inner(_, Some(itemType))))
          case (None, VArray(items)) => VArray(items.map(inner(_)))
          case (Some(TSchema(schemaName, fieldTypes)), VHash(fields)) =>
            val extra = fieldTypes.keySet.diff(fields.keySet)
            if (extra.nonEmpty) {
              throw ValueSerdeException(
                  s"invalid field(s) ${extra} in schema ${schemaName} value ${fields}"
              )
            }
            VHash(fieldTypes.collect {
              case (name, t) if fieldTypes.contains(name) =>
                name -> inner(fields(name), Some(t))
              case (name, t) if !Type.isOptional(t) =>
                throw ValueSerdeException(
                    s"missing non-optional member ${name} of schema ${schemaName}"
                )
            })
          case (_, VHash(fields)) =>
            VHash(fields.map { case (k, v) => k -> inner(v) })
          case (Some(TEnum(symbols)), s: VString) if symbols.contains(s.value) => s
          case (Some(TEnum(symbols)), other) =>
            throw ValueSerdeException(
                s"${other} is not one of the allowed symbols ${symbols.mkString(",")}"
            )
          case _ => innerValue
        }
      }
    }
    inner(value, t)
  }
}
