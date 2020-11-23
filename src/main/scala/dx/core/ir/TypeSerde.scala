package dx.core.ir

import dx.core.ir.Type._
import spray.json._
import dx.util.JsUtils

object TypeSerde {

  /**
    * Serialize a single type.
    * @param t the type to serialize
    * @param schemas schemas that t might reference
    * @return (serialized type, updated serialized schemas)
    */
  def serialize(
      t: Type,
      schemas: Map[String, JsValue] = Map.empty
  ): (JsValue, Map[String, JsValue]) = {
    t match {
      case TBoolean   => (JsString("Boolean"), schemas)
      case TInt       => (JsString("Int"), schemas)
      case TFloat     => (JsString("Float"), schemas)
      case TString    => (JsString("String"), schemas)
      case TFile      => (JsString("File"), schemas)
      case TDirectory => (JsString("Directory"), schemas)
      case THash      => (JsString("Hash"), schemas)
      case TArray(memberType, nonEmpty) =>
        val (typeJs, newSchemas) = serialize(memberType, schemas)
        (JsObject(
             Map(
                 "name" -> JsString("Array"),
                 "type" -> typeJs,
                 "nonEmpty" -> JsBoolean(nonEmpty)
             )
         ),
         newSchemas)
      case TSchema(name, _) if schemas.contains(name) =>
        (JsString(name), schemas)
      case TSchema(name, members) =>
        val (membersJs, newSchemas) =
          members.foldLeft((Map.empty[String, JsValue], Map.empty[String, JsValue])) {
            case ((membersAccu, aliasesAccu), (name, t)) =>
              val (typeJs, newSchemas) = serialize(t, aliasesAccu)
              (membersAccu + (name -> typeJs), newSchemas)
          }
        val schemaJs = JsObject(
            Map(
                "name" -> JsString(name),
                "members" -> JsObject(membersJs)
            )
        )
        (JsString(name), newSchemas + (name -> schemaJs))
      case TOptional(inner) =>
        serialize(inner, schemas) match {
          case (name: JsString, newSchemas) =>
            (JsObject(Map("name" -> name, "optional" -> JsBoolean(true))), newSchemas)
          case (JsObject(fields), newSchemas) =>
            (JsObject(fields + ("optional" -> JsBoolean(true))), newSchemas)
          case (other, _) =>
            throw new Exception(s"invalid inner type value ${other}")
        }
    }
  }

  /**
    * Serializes a mapping of variable names to Types.
    * @param types mapping of variable names to Types
    * @return a JsObject with two fields: 'types' and 'schemas'. Any TSchema in the
    *         input map are serialized to a JsString in the 'types' field and a
    *         corresponding entry in the 'schemas' field.
    */
  def serializeMap(
      types: Map[String, Type],
      jsSchema: Map[String, JsValue] = Map.empty,
      encodeDots: Boolean = true
  ): (Map[String, JsValue], Map[String, JsValue]) = {
    types.foldLeft((Map.empty[String, JsValue], jsSchema)) {
      case ((typeAccu, schemaAccu), (name, t)) =>
        val nameEncoded = if (encodeDots) {
          Parameter.encodeDots(name)
        } else {
          name
        }
        val (typeJs, newSchemas) = serialize(t, schemaAccu)
        (typeAccu + (nameEncoded -> typeJs), newSchemas)
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
            "schemas" -> JsObject(schemasJs)
        )
    )
  }

  /**
    * Serialize and create a specification for a single type.
    * @param t the type
    * @param schemas schemas that may be referenced by `t`
    * @return JsObject containing the specification
    */
  def serializeOne(t: Type, schemas: Map[String, JsValue] = Map.empty): JsObject = {
    val (typeJs, newSchemas) = serialize(t, schemas)
    JsObject("type" -> typeJs, "schemas" -> JsObject(newSchemas))
  }

  private def deserializeSchema(jsSchema: JsValue,
                                schemas: Map[String, TSchema],
                                jsSchemas: Map[String, JsValue]): Map[String, TSchema] = {
    jsSchema.asJsObject.getFields("name", "members") match {
      case Seq(JsString(name), JsObject(membersJs)) =>
        val (memberTypes, newSchemas) =
          membersJs.foldLeft((Map.empty[String, Type], schemas)) {
            case ((memberAccu, schemaAccu), (name, jsType)) =>
              val (t, newSchemas) = deserialize(jsType, schemaAccu, jsSchemas)
              (memberAccu + (name -> t), newSchemas)
          }
        newSchemas + (name -> TSchema(name, memberTypes))
      case _ =>
        throw new Exception(s"invalid schema ${jsSchema}")
    }
  }

  /**
    * Deserializes a serialized type value.
    * @param jsValue the serialized values
    * @param schemas schemas that the value may reference
    * @param jsSchemas serialized schemas that we only deserialize if
    *                  they are referenced
    * @return (type, new schema map)
    */
  def deserialize(
      jsValue: JsValue,
      schemas: Map[String, TSchema] = Map.empty,
      jsSchemas: Map[String, JsValue] = Map.empty
  ): (Type, Map[String, TSchema]) = {
    jsValue match {
      case JsString(name) if schemas.contains(name) =>
        (schemas(name), schemas)
      case JsString(name) if jsSchemas.contains(name) =>
        val newSchemas = deserializeSchema(jsSchemas(name), schemas, jsSchemas)
        (newSchemas(name), newSchemas)
      case JsString(name) =>
        (simpleFromString(name), schemas)
      case JsObject(fields) =>
        val (t, newSchemas) = fields("name") match {
          case JsString("Array") =>
            val (arrayType, newSchemas) = deserialize(fields("type"), schemas, jsSchemas)
            val nonEmpty = fields.get("nonEmpty").exists(JsUtils.getBoolean(_))
            (TArray(arrayType, nonEmpty), newSchemas)
          case JsString(name) if schemas.contains(name) =>
            (schemas(name), schemas)
          case JsString(name) if jsSchemas.contains(name) =>
            val newSchemas = deserializeSchema(jsSchemas(name), schemas, jsSchemas)
            (newSchemas(name), newSchemas)
          case JsString(name) =>
            (simpleFromString(name), schemas)
          case _ =>
            throw new Exception(s"invalid type field value ${jsValue}")
        }
        if (fields.get("optional").exists(JsUtils.getBoolean(_))) {
          (TOptional(t), newSchemas)
        } else {
          (t, newSchemas)
        }
      case _ =>
        throw new Exception(s"unexpected type value ${jsValue}")
    }
  }

  /**
    * Deserializes a map on parameter names to serialized values.
    * @param jsTypes types to deserialize
    * @param schemas schemas that may be referenced by the types
    * @param jsSchemas serialized schemas that we only deserialize if they are referenced
    * @param decodeDots whether to decode dots in parameter names
    * @return
    */
  def deserializeMap(
      jsTypes: Map[String, JsValue],
      schemas: Map[String, TSchema] = Map.empty,
      jsSchemas: Map[String, JsValue] = Map.empty,
      decodeDots: Boolean = true
  ): (Map[String, Type], Map[String, TSchema]) = {
    jsTypes.foldLeft((Map.empty[String, Type], schemas)) {
      case ((typeAccu, schemaAccu), (name, jsType)) =>
        val nameDecoded = if (decodeDots) {
          Parameter.decodeDots(name)
        } else {
          name
        }
        val (t, newSchemas) = deserialize(jsType, schemaAccu, jsSchemas)
        (typeAccu + (nameDecoded -> t), newSchemas)
    }
  }

  /**
    * Deserializes a parameter specification that was serialized using the
    * `serializeSpec` function.
    * @param jsValue the value to deserialize
    * @param schemas initial set of schemas (i.e. type aliases)
    * @param decodeDots whether to decode dots in variable names
    * @return mapping of variable names to deserialized Types
    */
  def deserializeSpec(jsValue: JsValue,
                      schemas: Map[String, TSchema] = Map.empty,
                      decodeDots: Boolean = true): Map[String, Type] = {
    val (jsTypes, jsSchemas) = jsValue match {
      case obj: JsObject if obj.fields.contains("types") =>
        obj.getFields("types", "schemas") match {
          case Seq(JsObject(jsTypes), JsObject(jsAliases)) =>
            (jsTypes, jsAliases)
          case Seq(JsObject(jsTypes)) =>
            (jsTypes, Map.empty[String, JsValue])
          case _ =>
            throw new Exception(s"invalid serialized types value ${jsValue}")
        }
      case JsObject(jsTypes) =>
        (jsTypes, Map.empty[String, JsValue])
      case _ =>
        throw new Exception(s"invalid serialized types value ${jsValue}")
    }
    val (types, _) = deserializeMap(jsTypes, schemas, jsSchemas, decodeDots)
    types
  }

  /**
    * Deserialize a single JsValue that was serialized using the `serializeOne` function.
    * @param jsValue the value to deserialize
    * @param schemas initial set of schemas (i.e. type aliases)
    * @return
    */
  def deserializeOne(jsValue: JsValue,
                     schemas: Map[String, TSchema] = Map.empty): (Type, Map[String, TSchema]) = {
    val (jsType, jsSchemas) = jsValue match {
      case obj: JsObject if obj.fields.contains("type") =>
        obj.getFields("type", "schemas") match {
          case Seq(jsType, JsObject(jsAliases)) => (jsType, jsAliases)
          case Seq(jsType)                      => (jsType, Map.empty[String, JsValue])
          case _ =>
            throw new Exception(s"invalid serialized type value ${jsValue}")
        }
      case _ => (jsValue, Map.empty[String, JsValue])
    }
    deserialize(jsType, schemas, jsSchemas)
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

  def fromNativeSpec(fields: Vector[JsValue]): Map[String, Type] = {
    fields.map {
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
    }.toMap
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
