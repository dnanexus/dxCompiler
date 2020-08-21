package dx.core.ir

sealed trait Value

object Value {
  // Primitive values.
  case class VInt(value: Long) extends Value
  case class VFloat(value: Double) extends Value
  case class VString(value: String) extends Value
  case class VBoolean(value: Boolean) extends Value
  case class VFile(value: String) extends Value

  // Represents a DNAnexus folder or a file-based representation of a directory
  // (e.g. a zip or tar archive file).
  case class VDirectory(value: String) extends Value

  /**
    * Represents the empty value for an optional field.
    */
  case object VNull extends Value

  /**
    * An array of values
    */
  case class VArray(value: Vector[Value]) extends Value

  /**
    * A JSON object.
    */
  case class VHash(value: Map[String, Value]) extends Value

  /**
    * A JSON object that represents a map with arbitrary key type.
    * @param value the map values
    */
  case class VMap(value: Map[Value, Value]) extends Value
}
