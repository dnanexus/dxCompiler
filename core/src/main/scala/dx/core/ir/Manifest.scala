package dx.core.ir

import dx.api.{DxApi, DxFile, DxFileDescCache}
import dx.core.ir.Type._
import dx.core.ir.Value.VFile
import dx.util.JsUtils
import spray.json._

import java.nio.file.Path
import scala.collection.immutable.{SortedMap, TreeMap}

object Manifest {
  val IdKey = "id"
  val ValuesKey = "values"
  val TypesKey = "types"
  val DefinitionsKey = "definitions"

  val FullFormatKeys = Set(IdKey, TypesKey, ValuesKey)
  val AllKeys = Set(IdKey, ValuesKey, TypesKey, DefinitionsKey)
}

case class Manifest(values: Map[String, Value],
                    types: Map[String, Type],
                    definitions: Map[String, Type] = Map.empty,
                    id: Option[String] = None) {
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
    val fields = TreeMap(
        Manifest.DefinitionsKey -> JsObject(jsDefinitions2),
        Manifest.TypesKey -> JsObject(jsTypes),
        Manifest.ValuesKey -> JsObject(jsValues)
    )
    val idField =
      id.map(i => TreeMap(Manifest.IdKey -> JsString(i))).getOrElse(TreeMap.empty[String, JsValue])
    JsObject(fields ++ idField)
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

  private def parseValues(jsValues: Map[String, JsValue],
                          types: Map[String, Type]): Map[String, Value] = {
    types.flatMap {
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
  }

  def parseUserManifest(jsValue: JsValue, expectedTypes: Map[String, Type]): Manifest = {
    val values = parseValues(jsValue.asJsObject.fields, expectedTypes)
    Manifest(values, expectedTypes, Type.collectSchemas(expectedTypes.values.toVector))
  }

  def parseFullManifest(jsValue: JsValue,
                        expectedTypes: Option[Map[String, Type]] = None): Manifest = {
    val fields = jsValue.asJsObject.fields
    val id = fields.get(Manifest.IdKey) match {
      case Some(JsString(id)) => Some(id)
      case None               => None
      case other              => throw new Exception(s"invalid manifest ID ${other}")
    }
    val jsValues = fields.get(Manifest.ValuesKey) match {
      case Some(JsObject(values)) => values
      case other                  => throw new Exception(s"invalid manifest values ${other}")
    }
    val jsDefinitions: Map[String, JsValue] = fields.get(Manifest.DefinitionsKey) match {
      case Some(JsObject(fields)) => fields
      case None                   => Map.empty[String, JsValue]
      case other                  => throw new Exception(s"invalid 'definitions' ${other}")
    }
    val (types, definitions) = fields.get(Manifest.TypesKey) match {
      case Some(JsObject(typeDefs)) =>
        val (types, definitions) = TypeSerde.deserializeMap(typeDefs, typeAliases, jsDefinitions)
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
        (types, definitions)
      case None if expectedTypes.isDefined =>
        (expectedTypes.get, Map.empty[String, Type])
      case other => throw new Exception(s"invalid 'types' ${other}")
    }
    val values = parseValues(jsValues, types)
    Manifest(values, types, definitions, id)
  }

  def parse(jsValue: JsValue, expectedTypes: Option[Map[String, Type]] = None): Manifest = {
    jsValue match {
      case JsObject(fields) if Manifest.FullFormatKeys.forall(fields.contains) =>
        parseFullManifest(jsValue, expectedTypes)
      case JsObject(fields)
          if fields.contains(Manifest.ValuesKey) && fields.keySet
            .diff(Manifest.AllKeys)
            .isEmpty && expectedTypes.isDefined =>
        parseFullManifest(jsValue, expectedTypes)
      case _ if expectedTypes.isDefined =>
        parseUserManifest(jsValue, expectedTypes.get)
      case _ =>
        throw new Exception(s"type information is required to parse manifest ${jsValue}")
    }
  }

  def parseFile(path: Path, expectedTypes: Option[Map[String, Type]] = None): Manifest = {
    parse(JsUtils.jsFromFile(path), expectedTypes)
  }
}
