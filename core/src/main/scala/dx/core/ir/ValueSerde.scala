package dx.core.ir

import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.CollectionUtils.IterableOnceExtensions
import spray.json._

import scala.collection.immutable.TreeSeqMap

object ValueSerde extends DefaultJsonProtocol {
  val WrappedValueKey = "wrapped___"

  case class ValueSerdeException(message: String) extends Exception(message)

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

  def wrapValue(jsValue: JsValue, mustNotBeWrapped: Boolean = false): JsValue = {
    jsValue match {
      case JsObject(fields) if !fields.contains(WrappedValueKey) =>
        JsObject(WrappedValueKey -> jsValue)
      case _ if mustNotBeWrapped =>
        throw ValueSerdeException(s"expected ${jsValue} to not be wrapped")
      case _ => jsValue
    }
  }

  def serializeWithType(
      value: Value,
      irType: Type,
      handler: Option[(Value, Type) => Either[Value, JsValue]] = None
  ): JsValue = {
    lazy val handlerWrapper = handler.map(h => (v: Value) => h(v, TMulti.Any))
    def inner(innerValue: Value, innerType: Type): JsValue = {
      val v = handler.map(_(innerValue, innerType)) match {
        case Some(Right(result))  => return result
        case Some(Left(newValue)) => newValue
        case None                 => innerValue
      }
      (innerType, v) match {
        case (_: TOptional, VNull)                    => JsNull
        case (TOptional(t), _)                        => inner(v, t)
        case (t: TMulti, _)                           => innerMulti(v, t)
        case (TBoolean, VBoolean(b))                  => JsBoolean(b)
        case (TInt, VInt(i))                          => JsNumber(i)
        case (TInt, VFloat(i)) if i.isValidInt        => JsNumber(i.intValue())
        case (TFloat, VFloat(f))                      => JsNumber(f)
        case (TFloat, VInt(i))                        => JsNumber(i.floatValue())
        case (TString, VString(s))                    => JsString(s)
        case (TFile | TString, VFile(path))           => JsString(path)
        case (TDirectory | TString, VDirectory(path)) => JsString(path)
        case (TArray(_, true), VArray(items)) if items.isEmpty =>
          throw ValueSerdeException(s"empty value for non-empty array type ${innerType}")
        case (TArray(itemType, _), VArray(items)) =>
          JsArray(items.map(inner(_, itemType)))
        case (TSchema(schemaName, fieldTypes), VHash(fields)) =>
          val extra = fieldTypes.keySet.diff(fields.keySet)
          if (extra.nonEmpty) {
            throw ValueSerdeException(
                s"invalid field(s) ${extra} in schema ${schemaName} value ${fields}"
            )
          }
          JsObject(fieldTypes.collect {
            case (name, t) if fields.contains(name) => name -> inner(fields(name), t)
            case (name, t) if !Type.isOptional(t) =>
              throw new Exception(s"missing non-optional member ${name} of schema ${schemaName}")
          })
        case (THash, VHash(fields)) =>
          JsObject(fields.view.mapValues(inner(_, TMulti.Any)).toMap)
        case (TEnum(symbols), VString(s)) if symbols.contains(s) => JsString(s)
        case _ =>
          throw ValueSerdeException(s"cannot serialize ${innerValue} as ${innerType}")
      }
    }
    def innerMulti(innerValue: Value, innerMultiType: TMulti): JsValue = {
      val jsValue = if (innerMultiType.bounds.isEmpty) {
        serialize(innerValue, handlerWrapper)
      } else {
        innerMultiType.bounds.iterator
          .map { t =>
            try {
              Some(inner(innerValue, t))
            } catch {
              case _: Throwable => None
            }
          }
          .collectFirst {
            case Some(t) => t
          }
          .getOrElse(
              throw new Exception(s"value ${value} not serializable to ${innerMultiType.bounds}")
          )
      }
      wrapValue(jsValue)
    }
    inner(value, irType)
  }

  def serializeMap(values: Map[String, Value]): Map[String, JsValue] = {
    values.map {
      case (name, value) => name -> serialize(value)
    }
  }

  def isWrappedValue(jsValue: JsValue): Boolean = {
    jsValue match {
      case JsObject(fields) if fields.size == 1 && fields.contains(WrappedValueKey) => true
      case _                                                                        => false
    }
  }

  def unwrapValue(jsValue: JsValue, mustBeWrapped: Boolean = false): JsValue = {
    jsValue match {
      case JsObject(fields) if fields.size == 1 && fields.contains(WrappedValueKey) =>
        fields(WrappedValueKey)
      case _ if mustBeWrapped =>
        throw ValueSerdeException(s"not a wrapped value ${jsValue}")
      case _ => jsValue
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
        case _ if isWrappedValue(v)               => inner(unwrapValue(v))
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
        case (TOptional(_), JsNull) => VNull
        case (TOptional(t), _)      => inner(v, t)
        case (any: TMulti, _) if any.bounds.isEmpty && isWrappedValue(v) =>
          inner(unwrapValue(v), any)
        case (TMulti(bounds), _) if bounds.isEmpty => deserialize(v)
        case (TMulti(bounds), _) if isWrappedValue(v) =>
          val unwrappedValue = unwrapValue(v)
          bounds.iterator
            .collectFirstDefined { t =>
              try {
                Some(inner(unwrappedValue, t))
              } catch {
                case _: ValueSerdeException => None
              }
            }
            .getOrElse(
                throw ValueSerdeException(
                    s"value ${unwrappedValue} does not match any of ${bounds}"
                )
            )
        case (TBoolean, JsBoolean(b))                     => VBoolean(b.booleanValue)
        case (TInt, JsNumber(value)) if value.isValidLong => VInt(value.toLongExact)
        case (TFloat, JsNumber(value))                    => VFloat(value.toDouble)
        case (TString, JsString(s))                       => VString(s)
        case (TFile, JsString(path))                      => VFile(path)
        case (TDirectory, JsString(path))                 => VDirectory(path)
        case (TArray(_, true), JsArray(items)) if items.isEmpty =>
          throw ValueSerdeException(s"Cannot convert empty array to non-empty type ${innerType}")
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
            throw ValueSerdeException(
                s"struct ${name} value has members that do not appear in the struct definition: ${extra}"
            )
          }
          val missingNonOptional = keys1.diff(keys2).map(key => key -> fieldTypes(key)).filterNot {
            case (_, TOptional(_)) => false
            case _                 => true
          }
          if (missingNonOptional.nonEmpty) {
            throw ValueSerdeException(
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
          throw ValueSerdeException(s"cannot deserialize value ${innerValue} as type ${innerType}")
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
