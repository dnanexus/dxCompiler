package dx.core.ir

import dx.core.ir.Type._
import dx.core.ir.Value._
import spray.json._

import scala.collection.immutable.TreeSeqMap

object ValueSerde extends DefaultJsonProtocol {

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
          case (Some(TArray(_, true)), VArray(Vector())) =>
            throw new Exception("empty array for non-empty array type")
          case (Some(TArray(itemType, _)), VArray(items)) =>
            VArray(items.map(inner(_, Some(itemType))))
          case (None, VArray(items)) => VArray(items.map(inner(_)))
          case (Some(TSchema(name, fieldTypes)), VHash(fields)) =>
            VHash(fields.map {
              case (k, _) if !fieldTypes.contains(k) =>
                throw new Exception(s"invalid member ${k} of schema ${name}")
              case (k, v) => k -> inner(v, Some(fieldTypes(k)))
            })
          case (_, VHash(fields)) => VHash(fields.map { case (k, v) => k -> inner(v) })
          case (Some(TEnum(symbols)), s: VString) if symbols.contains(s.value) =>
            s
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

  /**
    * Serializes a Value to JSON.
    * @param value the Value to serialize
    * @param handler an optional function to perform special handling of certain values.
    *                If Right(jsValue) is returned, then jsValue is the result of the
    *                transformation. If Left(newValue) is returned, then newValue is
    *                transformed according to the default rules.
    * @return
    */
  def serialize(value: Value, handler: Option[Value => Either[Value, JsValue]] = None): JsValue = {
    def inner(innerValue: Value): JsValue = {
      val v = handler.map(_(innerValue)) match {
        case Some(Right(result))  => return result
        case Some(Left(newValue)) => newValue
        case None                 => innerValue
      }
      v match {
        case VNull            => JsNull
        case VBoolean(b)      => JsBoolean(b)
        case VInt(i)          => JsNumber(i)
        case VFloat(f)        => JsNumber(f)
        case VString(s)       => JsString(s)
        case VFile(path)      => JsString(path)
        case VDirectory(path) => JsString(path)
        case VArray(items)    => JsArray(items.map(inner))
        case VHash(fields)    => JsObject(fields.view.mapValues(inner).toMap)
      }
    }
    inner(value)
  }

  def serializeMap(values: Map[String, Value]): Map[String, JsValue] = {
    values.map {
      case (name, value) => name -> serialize(value)
    }
  }

  /**
    * Deserializes a JsValue to a Value, in the absence of type information.
    * @param jsValue the JsValue
    * @param handler an optional function to perform special handling of certain values.
    *                If Right(value) is returned, then value is the result of the
    *                transformation. If Left(newJsValue) is returned, then newJsValue is
    *                transformed according to the default rules.
    * @return
    */
  def deserialize(jsValue: JsValue,
                  handler: Option[JsValue => Either[JsValue, Value]] = None): Value = {
    def inner(innerValue: JsValue): Value = {
      val v = handler.map(_(innerValue)) match {
        case Some(Right(result))    => return result
        case Some(Left(newJsValue)) => newJsValue
        case None                   => innerValue
      }
      v match {
        case JsNull                               => VNull
        case JsBoolean(b)                         => VBoolean(b.booleanValue)
        case JsNumber(value) if value.isValidLong => VInt(value.toLongExact)
        case JsNumber(value)                      => VFloat(value.toDouble)
        case JsString(s)                          => VString(s)
        case JsArray(items)                       => VArray(items.map(x => inner(x)))
        case JsObject(fields)                     => VHash(fields.view.mapValues(inner).to(TreeSeqMap))
      }
    }
    inner(jsValue)
  }

  /**
    * Deserializes a JsValue to a Value of the specified type.
    * @param jsValue the JsValue
    * @param t the Type
    * @param handler an optional function to perform special handling of certain values.
    *                If Right(value) is returned, then value is the result of the
    *                transformation. If Left(newJsValue) is returned, then newJsValue is
    *                transformed according to the default rules.
    * @return
    */
  def deserializeWithType(
      jsValue: JsValue,
      t: Type,
      handler: Option[(JsValue, Type) => Either[JsValue, Value]] = None
  ): Value = {
    def inner(innerValue: JsValue, innerType: Type): Value = {
      val v = handler.map(_(innerValue, innerType)) match {
        case Some(Right(result))    => return result
        case Some(Left(newJsValue)) => newJsValue
        case None                   => innerValue
      }
      (innerType, v) match {
        case (TOptional(_), JsNull)                       => VNull
        case (TOptional(t), _)                            => inner(v, t)
        case (TBoolean, JsBoolean(b))                     => VBoolean(b.booleanValue)
        case (TInt, JsNumber(value)) if value.isValidLong => VInt(value.toLongExact)
        case (TFloat, JsNumber(value))                    => VFloat(value.toDouble)
        case (TString, JsString(s))                       => VString(s)
        case (TFile, JsString(path))                      => VFile(path)
        case (TDirectory, JsString(path))                 => VDirectory(path)
        case (TArray(_, true), JsArray(items)) if items.isEmpty =>
          throw new Exception(s"Cannot convert empty array to non-empty type ${innerType}")
        case (TArray(t, _), JsArray(items)) =>
          VArray(items.map(x => inner(x, t)))
        case (TSchema(name, fieldTypes), JsObject(fields)) =>
          // ensure 1) fields keys are a subset of typeTypes keys, 2) fields
          // values are convertable to the corresponding types, and 3) any keys
          // in fieldTypes that do not appear in fields are optional
          val keys1 = fields.keySet
          val keys2 = fieldTypes.keySet
          val extra = keys2.diff(keys1)
          if (extra.nonEmpty) {
            throw new Exception(
                s"struct ${name} value has members that do not appear in the struct definition: ${extra}"
            )
          }
          val missingNonOptional = keys1.diff(keys2).map(key => key -> fieldTypes(key)).filterNot {
            case (_, TOptional(_)) => false
            case _                 => true
          }
          if (missingNonOptional.nonEmpty) {
            throw new Exception(
                s"struct ${name} value is missing non-optional members ${missingNonOptional}"
            )
          }
          VHash(fieldTypes.collect {
            case (name, t) if fields.contains(name) => name -> inner(fields(name), t)
          })
        case (THash, JsObject(fields)) =>
          VHash(
              fields
                .map {
                  case (key, value) => key -> deserialize(value)
                }
                .to(TreeSeqMap)
          )
        case (TEnum(symbols), JsString(s)) if symbols.contains(s) =>
          VString(s)
        case _ =>
          throw new Exception(s"cannot deserialize value ${innerValue} as type ${innerType}")
      }
    }
    inner(jsValue, t)
  }

  def deserializeMap(m: Map[String, JsValue]): Map[String, Value] = {
    m.map {
      case (k, v) => k -> deserialize(v)
    }
  }

  // support automatic conversion to/from JsValue
  implicit val valueFormat: RootJsonFormat[Value] = new RootJsonFormat[Value] {
    override def read(jsv: JsValue): Value = deserialize(jsv)
    override def write(value: Value): JsValue = serialize(value)
  }
}
