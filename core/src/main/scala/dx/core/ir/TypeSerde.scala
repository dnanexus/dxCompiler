package dx.core.ir

import dx.core.ir.Type._
import spray.json._
import dx.util.JsUtils

import scala.collection.immutable.{SeqMap, SortedMap, TreeSeqMap}

/**
  * Functions for serialization and deserialization of types. Note that SortedMap is
  * used for all object types in JSON, which ensures that fields have a consistent
  * (lexicographic) ordering. This means that object values that are serialized and
  * then deserialized may have their fields reordered.
  */
object TypeSerde {
  def serializeSchema(
      t: TSchema,
      typeDefs: Map[String, JsValue] = Map.empty
  ): (JsValue, SortedMap[String, JsValue]) = {
    val (fieldsJs, newTypeDefs) =
      t.fields.foldLeft((SortedMap.empty[String, JsValue], typeDefs.to(SortedMap))) {
        case ((fieldsAccu, typeDefAccu), (name, t)) =>
          val (typeJs, newTypeDefs) = serialize(t, typeDefAccu)
          (fieldsAccu + (name -> typeJs), newTypeDefs)
      }
    (JsObject(
         "type" -> JsString(t.name),
         "fields" -> JsObject(fieldsJs)
     ),
     newTypeDefs)
  }

  /**
    * Serialize a single type.
    * @param t the type to serialize
    * @param typeDefs type definitions that might be referenced by t
    * @return (serialized type, updated serialized type defs)
    */
  def serialize(
      t: Type,
      typeDefs: Map[String, JsValue] = Map.empty
  ): (JsValue, SortedMap[String, JsValue]) = {
    val newTypeDefs = t match {
      case schemaType: TSchema if !typeDefs.contains(schemaType.name) =>
        val (schemaJs, newTypeDefs) = serializeSchema(schemaType, typeDefs)
        newTypeDefs + (schemaType.name -> schemaJs)
      case _ => typeDefs.to(SortedMap)
    }
    t match {
      case TBoolean         => (JsString("Boolean"), newTypeDefs)
      case TInt             => (JsString("Int"), newTypeDefs)
      case TFloat           => (JsString("Float"), newTypeDefs)
      case TString          => (JsString("String"), newTypeDefs)
      case TFile            => (JsString("File"), newTypeDefs)
      case TDirectory       => (JsString("Directory"), newTypeDefs)
      case THash            => (JsString("Hash"), newTypeDefs)
      case TSchema(name, _) => (JsString(name), newTypeDefs)
      case TArray(memberType, nonEmpty) =>
        val (typeJs, updatedTypeDefs) = serialize(memberType, newTypeDefs)
        (JsObject(
             "type" -> JsString("Array"),
             "items" -> typeJs,
             "nonEmpty" -> JsBoolean(nonEmpty)
         ),
         updatedTypeDefs)
      case TEnum(symbols) =>
        (JsObject(
             "type" -> JsString("Enum"),
             "symbols" -> JsArray(symbols.map(JsString(_)))
         ),
         newTypeDefs)
      case TOptional(inner) =>
        serialize(inner, newTypeDefs) match {
          case (name: JsString, updatedTypeDefs) =>
            (JsObject(SortedMap("type" -> name, "optional" -> JsBoolean(true))), updatedTypeDefs)
          case (JsObject(fields), updatedTypeDefs) =>
            (JsObject(fields + ("optional" -> JsBoolean(true))), updatedTypeDefs)
          case (other, _) =>
            throw new Exception(s"invalid inner type value ${other}")
        }
    }
  }

  /**
    * Serializes a mapping of variable names to Types.
    * @param types mapping of variable names to Types
    * @return (parameter types, updated type definitions)
    */
  def serializeMap(
      types: Map[String, Type],
      jsTypeDefs: Map[String, JsValue] = Map.empty,
      encodeDots: Boolean = true
  ): (SortedMap[String, JsValue], SortedMap[String, JsValue]) = {
    types.foldLeft((SortedMap.empty[String, JsValue], jsTypeDefs.to(SortedMap))) {
      case ((typeAccu, typeDefAccu), (name, t)) =>
        val nameEncoded = if (encodeDots) {
          Parameter.encodeDots(name)
        } else {
          name
        }
        val (typeJs, newTypeDefs) = serialize(t, typeDefAccu)
        (typeAccu + (nameEncoded -> typeJs), newTypeDefs)
    }
  }

  /**
    * Serializes a parameter specification.
    * @param parameters parameter types
    * @param encodeDots whether to encode dots in input names
    * @return JsObject containing the input specification
    */
  def serializeSpec(parameters: Map[String, Type], encodeDots: Boolean = true): JsValue = {
    val (typesJs, schemasJs) = serializeMap(parameters, encodeDots = encodeDots)
    JsObject(
        Map(
            "types" -> JsObject(typesJs),
            "definitions" -> JsObject(schemasJs)
        )
    )
  }

  /**
    * Serialize and create a specification for a single type.
    * @param t the type
    * @param typeDefs type definitions that may be referenced by `t`
    * @return JsObject containing the specification
    */
  def serializeOne(t: Type, typeDefs: Map[String, JsValue] = Map.empty): JsObject = {
    val (typeJs, newTypeDefs) = serialize(t, typeDefs)
    JsObject("type" -> typeJs, "definitions" -> JsObject(newTypeDefs))
  }

  private def deserializeSchema(jsSchema: JsValue,
                                typeDefs: Map[String, Type],
                                jsTypeDefs: Map[String, JsValue],
                                name: Option[String] = None): Map[String, Type] = {
    val (schemaName, fieldsJs) = jsSchema.asJsObject.getFields("fields", "type") match {
      case Seq(JsObject(fieldsJs), JsString(name))   => (name, fieldsJs)
      case Seq(JsObject(fieldsJs)) if name.isDefined => (name.get, fieldsJs)
      case _ =>
        throw new Exception(s"invalid schema ${jsSchema}")
    }
    val (fieldTypes, newTypeDefs) =
      fieldsJs.foldLeft((Map.empty[String, Type], typeDefs)) {
        case ((fieldAccu, typeDefAccu), (name, jsType)) =>
          val (t, newTypeDefs) = deserialize(jsType, typeDefAccu, jsTypeDefs)
          (fieldAccu + (name -> t), newTypeDefs)
      }
    newTypeDefs + (schemaName -> TSchema(schemaName, fieldTypes.to(TreeSeqMap)))
  }

  def deserializeSchemas(
      jsSchemas: Map[String, JsValue],
      typeDefs: Map[String, Type] = Map.empty
  ): Map[String, Type] = {
    jsSchemas.values
      .foldLeft(typeDefs) {
        case (schemaAccu, jsSchema) => deserializeSchema(jsSchema, schemaAccu, jsSchemas)
      }
  }

  /**
    * Deserializes a serialized type value.
    * @param jsValue the serialized values
    * @param typeDefs type definitions that the may be referenced by jsValue
    * @param jsTypeDefs serialized type definitions that we only deserialize if
    *                  they are referenced
    * @return (type, updated type definition map)
    */
  def deserialize(
      jsValue: JsValue,
      typeDefs: Map[String, Type] = Map.empty,
      jsTypeDefs: Map[String, JsValue] = Map.empty
  ): (Type, Map[String, Type]) = {
    jsValue match {
      case JsString(name) if typeDefs.contains(name) =>
        (typeDefs(name), typeDefs)
      case JsString(name) if jsTypeDefs.contains(name) =>
        jsTypeDefs(name) match {
          case obj: JsObject if obj.fields.contains("type") =>
            deserialize(obj, typeDefs, jsTypeDefs)
          case obj: JsObject =>
            val newTypeDefs = deserializeSchema(obj, typeDefs, jsTypeDefs, Some(name))
            (newTypeDefs(name), newTypeDefs)
        }
      case JsString(name) =>
        (simpleFromString(name), typeDefs)
      case JsObject(fields) =>
        val (t, newTypeDefs) = fields("type") match {
          case JsString("Array") =>
            val (arrayType, newTypeDefs) = deserialize(fields("items"), typeDefs, jsTypeDefs)
            val nonEmpty = fields.get("nonEmpty").exists(JsUtils.getBoolean(_))
            (TArray(arrayType, nonEmpty), newTypeDefs)
          case JsString("Enum") =>
            val symbols = fields("symbols") match {
              case JsArray(values) =>
                values.map {
                  case JsString(s) => s
                  case other       => throw new Exception(s"Invalid enum symbol ${other}")
                }
              case other => throw new Exception(s"Invalid enum symbols ${other}")
            }
            (TEnum(symbols), typeDefs)
          case JsString(name) if typeDefs.contains(name) =>
            (typeDefs(name), typeDefs)
          case JsString(name) if jsTypeDefs.contains(name) =>
            val newTypeDefs = deserializeSchema(jsTypeDefs(name), typeDefs, jsTypeDefs)
            (newTypeDefs(name), newTypeDefs)
          case JsString(name) =>
            (simpleFromString(name), typeDefs)
          case _ =>
            throw new Exception(s"invalid type field value ${jsValue}")
        }
        if (fields.get("optional").exists(JsUtils.getBoolean(_))) {
          (TOptional(t), newTypeDefs)
        } else {
          (t, newTypeDefs)
        }
      case _ =>
        throw new Exception(s"unexpected type value ${jsValue}")
    }
  }

  /**
    * Deserializes a map on parameter names to serialized values.
    * @param jsTypes types to deserialize
    * @param typeDefs type definitions that may be referenced by the types
    * @param jsTypeDefs serialized type definitions that we only deserialize
    *                   if they are referenced
    * @param decodeDots whether to decode dots in parameter names
    * @return (parameter types, updated type definitions)
    */
  def deserializeMap(
      jsTypes: Map[String, JsValue],
      typeDefs: Map[String, Type] = Map.empty,
      jsTypeDefs: Map[String, JsValue] = Map.empty,
      decodeDots: Boolean = true
  ): (Map[String, Type], Map[String, Type]) = {
    jsTypes.foldLeft((Map.empty[String, Type], typeDefs)) {
      case ((typeAccu, typeDefAccu), (name, jsType)) =>
        val nameDecoded = if (decodeDots) {
          Parameter.decodeDots(name)
        } else {
          name
        }
        val (t, newTypeDefs) = deserialize(jsType, typeDefAccu, jsTypeDefs)
        (typeAccu + (nameDecoded -> t), newTypeDefs)
    }
  }

  /**
    * Deserializes a parameter specification that was serialized using the
    * `serializeSpec` function.
    * @param jsValue the value to deserialize
    * @param typeDefs initial set of schemas (i.e. type aliases)
    * @param decodeDots whether to decode dots in variable names
    * @return mapping of variable names to deserialized Types
    */
  def deserializeSpec(jsValue: JsValue,
                      typeDefs: Map[String, TSchema] = Map.empty,
                      decodeDots: Boolean = true): Map[String, Type] = {
    val (jsTypes, jsTypeDefs) = jsValue match {
      case obj: JsObject if obj.fields.contains("types") =>
        obj.getFields("types", "definitions") match {
          case Seq(JsObject(jsTypes), JsObject(jsDefinitions)) =>
            (jsTypes, jsDefinitions)
          case Seq(JsObject(jsTypes)) =>
            (jsTypes, Map.empty[String, JsValue])
          case _ =>
            throw new Exception(s"invalid serialized types or definitions in ${jsValue}")
        }
      case JsObject(jsTypes) =>
        (jsTypes, Map.empty[String, JsValue])
      case _ =>
        throw new Exception(s"invalid serialized spec ${jsValue}")
    }
    val (types, _) = deserializeMap(jsTypes, typeDefs, jsTypeDefs, decodeDots)
    types
  }

  /**
    * Deserialize a single JsValue that was serialized using the `serializeOne` function.
    * @param jsValue the value to deserialize
    * @param typeDefs initial set of type definitions
    * @return (deserialized type, updated type definitions)
    */
  def deserializeOne(jsValue: JsValue,
                     typeDefs: Map[String, Type] = Map.empty): (Type, Map[String, Type]) = {
    val (jsType, jsTypeDefs) = jsValue match {
      case obj: JsObject if obj.fields.contains("type") =>
        obj.getFields("type", "definitions") match {
          case Seq(jsType, JsObject(jsDefinitions)) => (jsType, jsDefinitions)
          case Seq(jsType)                          => (jsType, Map.empty[String, JsValue])
          case _ =>
            throw new Exception(s"invalid serialized type or definitions in ${jsValue}")
        }
      case _ => (jsValue, Map.empty[String, JsValue])
    }
    deserialize(jsType, typeDefs, jsTypeDefs)
  }

  private def toNativePrimitive(t: Type): String = {
    t match {
      case TBoolean => "boolean"
      case TInt     => "int"
      case TFloat   => "float"
      case TString  => "string"
      case TFile    => "file"
      // TODO: case TDirectory =>
      case _ => throw new Exception(s"not a primitive type")
    }
  }

  def toNative(t: Type): (String, Boolean) = {
    val (innerType, optional) = t match {
      case TOptional(innerType) => (innerType, true)
      case _                    => (t, false)
    }
    innerType match {
      case _ if Type.isNativePrimitive(innerType) =>
        (toNativePrimitive(innerType), optional)
      case TArray(memberType, nonEmpty) if Type.isNativePrimitive(memberType) =>
        // arrays of primitives translate to e.g. 'array:file' -
        val nativeInnerType = toNativePrimitive(memberType)
        (s"array:${nativeInnerType}", !nonEmpty || optional)
      case _ =>
        // everything else is a complex type represented as a hash
        ("hash", optional)
    }
  }

  private def fromNativeNonArray(cls: String): Type = {
    cls match {
      case "boolean" => TBoolean
      case "int"     => TInt
      case "float"   => TFloat
      case "string"  => TString
      case "file"    => TFile
      case "hash"    => THash
      case _         => throw new Exception(s"invalid native class ${cls}")
    }
  }

  def fromNative(cls: String, optional: Boolean): Type = {
    val t = if (cls.startsWith("array:")) {
      TArray(fromNativeNonArray(cls.drop(6)))
    } else {
      fromNativeNonArray(cls)
    }
    if (optional) {
      TOptional(t)
    } else {
      t
    }
  }

  def fromNativeSpec(fields: Vector[JsValue]): SeqMap[String, Type] = {
    fields
      .map {
        case JsObject(field) =>
          val name = field.get("name") match {
            case Some(JsString(name)) => name
            case _ =>
              throw new Exception(s"invalid or missing name for field ${fields}")
          }
          val cls = field.get("class") match {
            case Some(JsString(cls)) => cls
            case None                => "string"
            case other =>
              throw new Exception(s"invalid native 'class' value ${other}")
          }
          val optional = field.get("optional") match {
            case Some(JsBoolean(optional)) => optional
            case None                      => false
            case other =>
              throw new Exception(s"invalid native 'optional' value ${other}")
          }
          name -> fromNative(cls, optional)
        case other =>
          throw new Exception(s"invalid native input/output spec field ${other}")
      }
      .to(TreeSeqMap)
  }

  // Get a human readable type name
  // Int ->   "Int"
  // Array[Int] -> "Array[Int]"
  def toString(t: Type): String = {
    t match {
      case TBoolean         => "Boolean"
      case TInt             => "Int"
      case TFloat           => "Float"
      case TString          => "String"
      case TFile            => "File"
      case TDirectory       => "Directory"
      case THash            => "Hash"
      case TSchema(name, _) => name
      case TArray(memberType, _) =>
        s"Array[${toString(memberType)}]"
      case TEnum(allowedValues) =>
        s"Enum{${allowedValues.mkString(",")}}"
      case TOptional(TOptional(_)) =>
        throw new Exception(s"nested optional type ${t}")
      case TOptional(inner) =>
        s"${toString(inner)}?"
    }
  }

  case class UnknownTypeException(message: String) extends Exception(message)

  /**
    * Convert a String to a simple (non-compound) type, i.e. TArray and TMap
    * are not supported.
    * @param s type string
    * @return
    */
  def simpleFromString(s: String): Type = {
    s match {
      case "Boolean"   => TBoolean
      case "Int"       => TInt
      case "Float"     => TFloat
      case "String"    => TString
      case "File"      => TFile
      case "Directory" => TDirectory
      case "Hash"      => THash
      case _ if s.endsWith("?") =>
        simpleFromString(s.dropRight(1)) match {
          case TOptional(_) =>
            throw new Exception(s"nested optional type ${s}")
          case inner =>
            TOptional(inner)
        }
      case s if s.contains("[") =>
        throw new Exception(s"type ${s} is not primitive")
      case _ =>
        throw UnknownTypeException(s"Unknown type ${s}")
    }
  }
}
