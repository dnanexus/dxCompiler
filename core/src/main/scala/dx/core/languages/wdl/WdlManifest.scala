package dx.core.languages.wdl

import dx.api.{DxApi, DxFile, DxFileDescCache}
import dx.util.JsUtils
import spray.json._
import wdlTools.eval.{WdlValueSerde, WdlValues}
import wdlTools.types.WdlTypes.T_Struct
import wdlTools.types.{TypeUtils, WdlTypeSerde, WdlTypes}

import java.nio.file.Path
import scala.collection.immutable.SortedMap

case class WdlManifest(values: Map[String, WdlValues.V],
                       types: Map[String, WdlTypes.T],
                       schemas: Map[String, WdlTypes.T]) {
  values.keySet.diff(types.keySet) match {
    case d if d.nonEmpty =>
      throw new Exception(s"missing type definition(s) ${d.mkString(",")}")
    case _ => ()
  }

  def toJson: JsObject = {
    val (jsTypes, allSchemas, jsSchemas) =
      types.foldLeft(SortedMap.empty[String, JsValue], schemas, SortedMap.empty[String, JsValue]) {
        case ((typeAccu, schemaAccu, schemaJsAccu), (name, t: WdlTypes.T_Struct))
            if schemas.contains(t.name) =>
          (typeAccu + (name -> JsString(t.name)), schemaAccu, schemaJsAccu)
        case ((typeAccu, schemaAccu, schemaJsAccu), (name, t: WdlTypes.T_Struct)) =>
          (typeAccu + (name -> JsString(t.name)), schemaAccu + (t.name -> t), schemaJsAccu)
        case ((typeAccu, schemaAccu, schemaJsAccu), (name, wdlType)) =>
          val (jsType, newSchemasJs) = WdlTypeSerde.serializeType(wdlType, schemaJsAccu)
          (typeAccu + (name -> jsType), schemaAccu, newSchemasJs)
      }
    val jsSchemas2 = allSchemas.foldLeft(jsSchemas) {
      case (accu, (name, _)) if accu.contains(name) => accu
      case (accu, (name, schema: WdlTypes.T_Struct)) =>
        val (schemaJs, newJsSchemas) =
          WdlTypeSerde.serializeStruct(schema, accu, withName = false)
        newJsSchemas + (name -> schemaJs)
    }
    val jsValues = WdlValueSerde.serializeMap(values)
    JsObject(
        "schemas" -> JsObject(jsSchemas2),
        "types" -> JsObject(jsTypes),
        "values" -> JsObject(jsValues)
    )
  }
}

case class WdlManifestParser(typeAliases: Map[String, WdlTypes.T] = Map.empty,
                             dxFileDescCache: DxFileDescCache = DxFileDescCache.empty,
                             dxApi: DxApi = DxApi.get) {
  private def fileHandler(value: JsValue, wdlType: WdlTypes.T): Option[WdlValues.V] = {
    (TypeUtils.unwrapOptional(wdlType), value) match {
      case (WdlTypes.T_File, fileObj: JsObject) if DxFile.isLinkJson(value) =>
        // Convert the dx link to a URI string. We can later decide if we want to download it or not.
        // Use the cache value if there is one to save the API call.
        val dxFile = dxFileDescCache.updateFileFromCache(DxFile.fromJson(dxApi, fileObj))
        Some(WdlValues.V_File(dxFile.asUri))
      case _ => None
    }
  }

  def parse(jsValue: JsValue,
            expectedTypes: Option[Map[String, WdlTypes.T]] = None): WdlManifest = {
    val (values, types, schemas) = jsValue match {
      case JsObject(fields) if fields.contains("types") && fields.contains("values") =>
        // extended manifest format
        val jsValues = fields("values") match {
          case JsObject(values) => values
          case other            => throw new Exception(s"Invalid manifest values ${other}")
        }
        val jsSchemas: Map[String, JsValue] = fields.get("schemas") match {
          case Some(JsObject(fields)) => fields
          case None                   => Map.empty
          case other                  => throw new Exception(s"invalid 'schemas' ${other}")
        }
        val (types, schemas) = fields("types") match {
          case JsObject(fields) => WdlTypeSerde.deserializeTypes(fields, typeAliases, jsSchemas)
          case other            => throw new Exception(s"invalid 'types' ${other}")
        }
        expectedTypes.foreach {
          case expected if expected.keySet != types.keySet =>
            throw new Exception(
                s"mismatch between actual and expected types: ${types} != ${expected}"
            )
          case expected =>
            expected.foreach {
              case (name, wdlType) if wdlType != types(name) =>
                throw new Exception(
                    s"mismatch between '${name}' actual and expected type definition: ${wdlType} != ${types(name)}"
                )
            }
        }
        val values = types.flatMap {
          case (name, wdlType) if jsValues.contains(name) =>
            Some(
                name -> WdlValueSerde
                  .deserializeWithType(jsValues(name), wdlType, name, Some(fileHandler))
            )
          case (_, _: WdlTypes.T_Optional) =>
            None
          case wdlType =>
            throw new Exception(s"missing value for non-optional type ${wdlType}")
        }
        (values, types, schemas.collect {
          case (name, struct: T_Struct) => name -> struct
        })
      case JsObject(jsValues) if expectedTypes.isDefined =>
        // simple manifest format
        val values = expectedTypes.get.flatMap {
          case (name, wdlType) if jsValues.contains(name) =>
            Some(
                name -> WdlValueSerde
                  .deserializeWithType(jsValues(name), wdlType, name, Some(fileHandler))
            )
          case (_, _: WdlTypes.T_Optional) =>
            None
          case wdlType =>
            throw new Exception(s"missing value for non-optional type ${wdlType}")
        }
        val schemas = TypeUtils.collectStructs(expectedTypes.get.values.toVector)
        (values, expectedTypes.get, schemas)
      case _ =>
        throw new Exception("expectedTypes must be specified for simple manifest")
    }
    WdlManifest(values, types, schemas)
  }

  def parseFile(path: Path, expectedTypes: Option[Map[String, WdlTypes.T]] = None): WdlManifest = {
    parse(JsUtils.jsFromFile(path), expectedTypes)
  }
}
