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
  val TypeKey = "type"
  val ItemsKey = "items"
  val NonEmptyKey = "nonEmpty"
  val FieldsKey = "fields"
  val SymbolsKey = "symbols"
  val ChoicesKey = "choices"
  val OptionalKey = "optional"
  val TypesKey = "types"
  val DefinitionsKey = "definitions"

  val ArrayTypeName = "Array"
  val EnumTypeName = "Enum"
  val MultiTypeName = "Multi"

  case class TypeSerdeException(message: String) extends Exception(message)

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
         TypeKey -> JsString(t.name),
         FieldsKey -> JsObject(fieldsJs)
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
             TypeKey -> JsString(ArrayTypeName),
             ItemsKey -> typeJs,
             NonEmptyKey -> JsBoolean(nonEmpty)
         ),
         updatedTypeDefs)
      case TEnum(symbols) =>
        (JsObject(
             TypeKey -> JsString(EnumTypeName),
             SymbolsKey -> JsArray(symbols.map(JsString(_)))
         ),
         newTypeDefs)
      case TOptional(inner) =>
        serialize(inner, newTypeDefs) match {
          case (name: JsString, updatedTypeDefs) =>
            (JsObject(TypeKey -> name, OptionalKey -> JsBoolean(true)), updatedTypeDefs)
          case (JsObject(fields), updatedTypeDefs) =>
            (JsObject(fields + (OptionalKey -> JsBoolean(true))), updatedTypeDefs)
          case (other, _) =>
            throw TypeSerdeException(s"invalid inner type value ${other}")
        }
      case TMulti(types) =>
        val (serializedTypes, updatedTypeDefs) =
          types.foldLeft(Vector.empty[JsValue], newTypeDefs) {
            case ((serializedTypeAccu, typeDefAccu), t) =>
              val (serializedType, updatedTypeDefs) = serialize(t, typeDefAccu)
              (serializedTypeAccu :+ serializedType, updatedTypeDefs)
          }
        (JsObject(TypeKey -> JsString(MultiTypeName), ChoicesKey -> JsArray(serializedTypes)),
         updatedTypeDefs)
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
          Parameter.encodeName(name)
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
            TypesKey -> JsObject(typesJs),
            DefinitionsKey -> JsObject(schemasJs)
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
    JsObject(TypeKey -> typeJs, DefinitionsKey -> JsObject(newTypeDefs))
  }

  private def deserializeSchema(jsSchema: JsValue,
                                typeDefs: Map[String, Type],
                                jsTypeDefs: Map[String, JsValue],
                                name: Option[String] = None): Map[String, Type] = {
    val (schemaName, fieldsJs) = jsSchema.asJsObject.getFields(FieldsKey, TypeKey) match {
      case Seq(JsObject(fieldsJs), JsString(name))   => (name, fieldsJs)
      case Seq(JsObject(fieldsJs)) if name.isDefined => (name.get, fieldsJs)
      case _ =>
        throw TypeSerdeException(s"invalid schema ${jsSchema}")
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
          case obj: JsObject if obj.fields.contains(TypeKey) =>
            deserialize(obj, typeDefs, jsTypeDefs)
          case obj: JsObject =>
            val newTypeDefs = deserializeSchema(obj, typeDefs, jsTypeDefs, Some(name))
            (newTypeDefs(name), newTypeDefs)
        }
      case JsString(name) =>
        (simpleFromString(name), typeDefs)
      case JsObject(fields) =>
        val (t, newTypeDefs) = fields(TypeKey) match {
          case JsString(ArrayTypeName) =>
            val (arrayType, newTypeDefs) = deserialize(fields(ItemsKey), typeDefs, jsTypeDefs)
            val nonEmpty = fields.get(NonEmptyKey).exists(JsUtils.getBoolean(_))
            (TArray(arrayType, nonEmpty), newTypeDefs)
          case JsString(EnumTypeName) =>
            val symbols = fields(SymbolsKey) match {
              case JsArray(values) =>
                values.map {
                  case JsString(s) => s
                  case other       => throw TypeSerdeException(s"invalid enum symbol ${other}")
                }
              case other => throw TypeSerdeException(s"invalid enum symbols ${other}")
            }
            (TEnum(symbols), typeDefs)
          case JsString(MultiTypeName) =>
            val (choices, newTypeDefs) = fields.get(ChoicesKey) match {
              case Some(JsArray(choices)) =>
                choices.foldLeft(Vector.empty[Type], typeDefs) {
                  case ((typeAccu, typeDefAccu), jsValue) =>
                    val (t, newTypeDefs) = deserialize(jsValue, typeDefAccu, jsTypeDefs)
                    (typeAccu :+ t, newTypeDefs)
                }
              case Some(JsNull) | None =>
                (Vector.empty[Type], typeDefs)
              case other =>
                throw TypeSerdeException(s"invalid multi-type array ${other}")
            }
            (TMulti(choices), newTypeDefs)
          case JsString(name) if typeDefs.contains(name) =>
            (typeDefs(name), typeDefs)
          case JsString(name) if jsTypeDefs.contains(name) =>
            val newTypeDefs = deserializeSchema(jsTypeDefs(name), typeDefs, jsTypeDefs)
            (newTypeDefs(name), newTypeDefs)
          case JsString(name) =>
            (simpleFromString(name), typeDefs)
          case _ =>
            throw TypeSerdeException(s"invalid type field value ${jsValue}")
        }
        if (fields.get(OptionalKey).exists(JsUtils.getBoolean(_))) {
          (TOptional(t), newTypeDefs)
        } else {
          (t, newTypeDefs)
        }
      case _ =>
        throw TypeSerdeException(s"unexpected type value ${jsValue}")
    }
  }

  /**
    * Deserializes a map on parameter names to serialized values.
    * @param jsTypes types to deserialize
    * @param typeDefs type definitions that may be referenced by the types
    * @param jsTypeDefs serialized type definitions that we only deserialize
    *                   if they are referenced
    * @param decodeNames whether to decode dots in parameter names
    * @return (parameter types, updated type definitions)
    */
  def deserializeMap(
      jsTypes: Map[String, JsValue],
      typeDefs: Map[String, Type] = Map.empty,
      jsTypeDefs: Map[String, JsValue] = Map.empty,
      decodeNames: Boolean = true
  ): (Map[String, Type], Map[String, Type]) = {
    jsTypes.foldLeft((Map.empty[String, Type], typeDefs)) {
      case ((typeAccu, typeDefAccu), (name, jsType)) =>
        val nameDecoded = if (decodeNames) {
          Parameter.decodeName(name)
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
    * @param decodeNames whether to decode dots in variable names
    * @return mapping of variable names to deserialized Types
    */
  def deserializeSpec(jsValue: JsValue,
                      typeDefs: Map[String, TSchema] = Map.empty,
                      decodeNames: Boolean = true): Map[String, Type] = {
    val (jsTypes, jsTypeDefs) = jsValue match {
      case obj: JsObject if obj.fields.contains(TypesKey) =>
        obj.getFields(TypesKey, DefinitionsKey) match {
          case Seq(JsObject(jsTypes), JsObject(jsDefinitions)) =>
            (jsTypes, jsDefinitions)
          case Seq(JsObject(jsTypes)) =>
            (jsTypes, Map.empty[String, JsValue])
          case _ =>
            throw TypeSerdeException(s"invalid serialized types or definitions in ${jsValue}")
        }
      case JsObject(jsTypes) =>
        (jsTypes, Map.empty[String, JsValue])
      case _ =>
        throw TypeSerdeException(s"invalid serialized spec ${jsValue}")
    }
    val (types, _) = deserializeMap(jsTypes, typeDefs, jsTypeDefs, decodeNames)
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
      case obj: JsObject if obj.fields.contains(TypeKey) =>
        obj.getFields(TypeKey, DefinitionsKey) match {
          case Seq(jsType, JsObject(jsDefinitions)) => (jsType, jsDefinitions)
          case Seq(jsType)                          => (jsType, Map.empty[String, JsValue])
          case _ =>
            throw TypeSerdeException(s"invalid serialized type or definitions in ${jsValue}")
        }
      case _ => (jsValue, Map.empty[String, JsValue])
    }
    deserialize(jsType, typeDefs, jsTypeDefs)
  }

  private def toNativePrimitive(t: Type, pathsAreNative: Boolean): String = {
    t match {
      case TBoolean                     => "boolean"
      case TInt                         => "int"
      case TFloat                       => "float"
      case TString                      => "string"
      case TFile if pathsAreNative      => "file"
      case TDirectory if pathsAreNative => "file"
      case TFile | TDirectory           => "hash"
      case _                            => throw TypeSerdeException(s"not a primitive type")
    }
  }

  def toNative(t: Type, pathsAreNative: Boolean = true): (String, Boolean) = {
    val (innerType, optional) = t match {
      case TOptional(innerType) => (innerType, true)
      case _                    => (t, false)
    }
    innerType match {
      case _ if Type.isNativePrimitive(innerType, pathsAreNative) =>
        (toNativePrimitive(innerType, pathsAreNative), optional)
      case TArray(memberType, nonEmpty) if Type.isNativePrimitive(memberType, pathsAreNative) =>
        // arrays of primitives translate to e.g. 'array:file' -
        val nativeInnerType = toNativePrimitive(memberType, pathsAreNative)
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
      case _         => throw TypeSerdeException(s"invalid native class ${cls}")
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
              throw TypeSerdeException(s"invalid or missing name for field ${fields}")
          }
          val cls = field.get("class") match {
            case Some(JsString(cls)) => cls
            case None                => "string"
            case other =>
              throw TypeSerdeException(s"invalid native 'class' value ${other}")
          }
          val optional = field.get("optional") match {
            case Some(JsBoolean(optional)) => optional
            case None                      => false
            case other =>
              throw TypeSerdeException(s"invalid native 'optional' value ${other}")
          }
          name -> fromNative(cls, optional)
        case other =>
          throw TypeSerdeException(s"invalid native input/output spec field ${other}")
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
        s"${ArrayTypeName}[${toString(memberType)}]"
      case TEnum(symbols) =>
        s"${EnumTypeName}{${symbols.mkString(",")}}"
      case TMulti(Vector()) => MultiTypeName
      case TMulti(choices) =>
        s"${MultiTypeName}{${choices.map(toString).mkString(",")}}"
      case TOptional(TOptional(_)) =>
        throw TypeSerdeException(s"nested optional type ${t}")
      case TOptional(inner) =>
        s"${toString(inner)}?"
    }
  }

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
            throw TypeSerdeException(s"nested optional type ${s}")
          case inner =>
            TOptional(inner)
        }
      case s if s.contains("[") =>
        throw TypeSerdeException(s"type ${s} is not primitive")
      case _ =>
        throw TypeSerdeException(s"Unknown type ${s}")
    }
  }
}
