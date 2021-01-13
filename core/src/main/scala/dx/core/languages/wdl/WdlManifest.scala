package dx.core.languages.wdl

import dx.api.{DxApi, DxFile, DxFileDescCache}
import dx.util.JsUtils
import spray.json._
import wdlTools.eval.{WdlValueSerde, WdlValues}
import wdlTools.types.{TypeUtils, WdlTypeSerde, WdlTypes}

import java.nio.file.Path

case class WdlManifest(values: Map[String, (WdlTypes.T, WdlValues.V)],
                       typeAliases: Map[String, WdlTypes.T]) {
  def serialize: JsObject = {}
}

case class WdlManifestParser(dxFileDescCache: DxFileDescCache, dxApi: DxApi = DxApi.get) {
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
    val (values, typeAliases) = jsValue match {
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
        val (wdlTypes, typeAliases) = fields("types") match {
          case JsObject(fields) => WdlTypeSerde.deserializeTypes(fields, jsSchemas = jsSchemas)
          case other            => throw new Exception(s"invalid 'types' ${other}")
        }
        expectedTypes.foreach {
          case expected if expected.keySet != wdlTypes.keySet =>
            throw new Exception(
                s"mismatch between actual and expected types: ${wdlTypes} != ${expected}"
            )
          case expected =>
            expected.foreach {
              case (name, wdlType) if wdlType != wdlTypes(name) =>
                throw new Exception(
                    s"mismatch between '${name}' actual and expected type definition: ${wdlType} != ${wdlTypes(name)}"
                )
            }
        }
        val values = wdlTypes.map {
          case (name, wdlType) if jsValues.contains(name) =>
            name -> (wdlType, WdlValueSerde.deserializeWithType(jsValues(name),
                                                                wdlType,
                                                                name,
                                                                Some(fileHandler)))
          case (name, wdlType: WdlTypes.T_Optional) =>
            name -> (wdlType, WdlValues.V_Null)
          case wdlType =>
            throw new Exception(s"missing value for non-optional type ${wdlType}")
        }
        (values, typeAliases)
      case JsObject(jsValues) if expectedTypes.isDefined =>
        // simple manifest format
        val values = expectedTypes.get.map {
          case (name, wdlType) if jsValues.contains(name) =>
            name -> (wdlType, WdlValueSerde.deserializeWithType(jsValues(name),
                                                                wdlType,
                                                                name,
                                                                Some(fileHandler)))
          case (name, wdlType: WdlTypes.T_Optional) =>
            name -> (wdlType, WdlValues.V_Null)
          case wdlType =>
            throw new Exception(s"missing value for non-optional type ${wdlType}")
        }
        (values, Map.empty[String, WdlTypes.T])
      case _ =>
        throw new Exception("expectedTypes must be specified for simple manifest")
    }
    WdlManifest(values, typeAliases)
  }

  def parseFile(path: Path, expectedTypes: Option[Map[String, WdlTypes.T]] = None): WdlManifest = {
    parse(JsUtils.jsFromFile(path), expectedTypes)
  }
}
