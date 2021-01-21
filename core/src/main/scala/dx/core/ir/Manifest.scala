package dx.core.ir

import dx.api.{DxApi, DxFile, DxFileDescCache}
import dx.core.ir.Type._
import dx.core.ir.Value.VFile
import dx.util.JsUtils
import spray.json._

import java.nio.file.Path
import scala.collection.immutable.SortedMap

case class Manifest(id: String,
                    values: Map[String, Value],
                    types: Map[String, Type],
                    definitions: Map[String, Type] = Map.empty) {
  values.keySet.diff(types.keySet) match {
    case d if d.nonEmpty =>
      throw new Exception(s"missing type definition(s) ${d.mkString(",")}")
    case _ => ()
  }

  def toJson: JsObject = {
    val (jsTypes, allDefinitions, jsDefinitions) =
      types.foldLeft(SortedMap.empty[String, JsValue],
                     definitions,
                     SortedMap.empty[String, JsValue]) {
        case ((typeAccu, defAccu, defJsAccu), (name, t: TSchema)) if defAccu.contains(t.name) =>
          (typeAccu + (name -> JsString(t.name)), defAccu, defJsAccu)
        case ((typeAccu, defAccu, defJsAccu), (name, t: TSchema)) =>
          (typeAccu + (name -> JsString(t.name)), defAccu + (t.name -> t), defJsAccu)
        case ((typeAccu, defAccu, defJsAccu), (name, t)) =>
          val (jsType, newDefJs) = TypeSerde.serialize(t, defJsAccu)
          (typeAccu + (name -> jsType), defAccu, newDefJs)
      }
    val jsDefinitions2 = allDefinitions.foldLeft(jsDefinitions) {
      case (accu, (name, _)) if accu.contains(name) => accu
      case (accu, (name, schema: TSchema)) =>
        val (schemaJs, newJsDefinitions) = TypeSerde.serializeSchema(schema, accu)
        newJsDefinitions + (name -> schemaJs)
    }
    val jsValues = ValueSerde.serializeMap(values)
    JsObject(
        "id" -> JsString(id),
        "definitions" -> JsObject(jsDefinitions2),
        "types" -> JsObject(jsTypes),
        "values" -> JsObject(jsValues)
    )
  }
}

case class ManifestParser(typeAliases: Map[String, Type] = Map.empty,
                          dxFileDescCache: DxFileDescCache = DxFileDescCache.empty,
                          dxApi: DxApi = DxApi.get) {
  private def fileHandler(jsValue: JsValue, t: Type): Either[JsValue, Value] = {
    (Type.unwrapOptional(t), jsValue) match {
      case (TFile, fileObj: JsObject) if DxFile.isLinkJson(jsValue) =>
        // Convert the dx link to a URI string. We can later decide if we want to download it or not.
        // Use the cache value if there is one to save the API call.
        val dxFile = dxFileDescCache.updateFileFromCache(DxFile.fromJson(dxApi, fileObj))
        Right(VFile(dxFile.asUri))
      case _ => Left(jsValue)
    }
  }

  def parse(jsValue: JsValue, expectedTypes: Option[Map[String, Type]] = None): Manifest = {
    val fields = jsValue.asJsObject.fields
    val id = fields.get("id") match {
      case Some(JsString(id)) => id
      case other              => throw new Exception(s"invalid manifest ID ${other}")
    }
    val jsValues = fields.get("values") match {
      case Some(JsObject(values)) => values
      case other                  => throw new Exception(s"invalid manifest values ${other}")
    }
    val jsDefinitions: Map[String, JsValue] = fields.get("definitions") match {
      case Some(JsObject(fields)) => fields
      case None                   => Map.empty
      case other                  => throw new Exception(s"invalid 'definitions' ${other}")
    }
    val (types, definitions) = fields.get("types") match {
      case Some(JsObject(typeDefs)) =>
        TypeSerde.deserializeMap(typeDefs, typeAliases, jsDefinitions)
      case other => throw new Exception(s"invalid 'types' ${other}")
    }
    expectedTypes.foreach {
      case expected if expected.keySet != types.keySet =>
        throw new Exception(
            s"mismatch between actual and expected types: ${types} != ${expected}"
        )
      case expected =>
        expected.foreach {
          case (name, t) if t != types(name) =>
            throw new Exception(
                s"mismatch between '${name}' actual and expected type definition: ${t} != ${types(name)}"
            )
        }
    }
    val values = types.flatMap {
      case (name, t) if jsValues.contains(name) =>
        Some(
            name -> ValueSerde
              .deserializeWithType(jsValues(name), t, Some(fileHandler))
        )
      case (_, _: TOptional) =>
        None
      case t =>
        throw new Exception(s"missing value for non-optional type ${t}")
    }
    Manifest(id, values, types, definitions)
  }

  def parseFile(path: Path, expectedTypes: Option[Map[String, Type]] = None): Manifest = {
    parse(JsUtils.jsFromFile(path), expectedTypes)
  }
}
