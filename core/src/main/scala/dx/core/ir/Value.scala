package dx.core.ir

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
}
