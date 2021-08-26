package dx.core.ir

import dx.api.{DxApi, DxFileDescCache}
import dx.core.ir.Type._
import dx.util.JsUtils
import spray.json._

import java.nio.file.Path
import scala.collection.immutable.{SortedMap, TreeMap}

object Manifest {
  val IdKey = "id"
  val ValuesKey = "values"
  val TypesKey = "types"
  val DefinitionsKey = "definitions"

  val AllKeys = Set(IdKey, ValuesKey, TypesKey, DefinitionsKey)

  def parseFullManifest(jsValue: JsValue,
                        irTypes: Option[Map[String, Type]] = None,
                        typeAliases: Map[String, Type] = Map.empty): Manifest = {
    val fields = jsValue.asJsObject.fields
    val id = fields.get(Manifest.IdKey) match {
      case Some(JsString(id)) => Some(id)
      case None               => None
      case other              => throw new Exception(s"invalid manifest ID ${other}")
    }
    val jsValues: Map[String, JsValue] = fields.get(Manifest.ValuesKey) match {
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
        irTypes.foreach {
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
        (Some(types), Some(definitions))
      case None =>
        (irTypes, None)
      case _ =>
        throw new Exception(s"invalid manifest ${jsValue}")
    }
    Manifest(jsValues, types, definitions, id)
  }

  def parse(jsValue: JsValue,
            types: Option[Map[String, Type]] = None,
            typeAliases: Map[String, Type] = Map.empty): Manifest = {
    jsValue match {
      case JsObject(fields)
          if fields.contains(Manifest.ValuesKey) &&
            fields.keySet.diff(Manifest.AllKeys).isEmpty =>
        parseFullManifest(jsValue, types, typeAliases)
      case _ =>
        Manifest(
            jsValue.asJsObject.fields,
            types,
            types.map(t => Type.collectSchemas(t.values.toVector))
        )
    }
  }

  def parseFile(path: Path,
                types: Option[Map[String, Type]] = None,
                typeAliases: Map[String, Type] = Map.empty): Manifest = {
    parse(JsUtils.jsFromFile(path), types, typeAliases)
  }
}

case class Manifest(jsValues: Map[String, JsValue],
                    types: Option[Map[String, Type]] = None,
                    definitions: Option[Map[String, Type]] = None,
                    id: Option[String] = None) {
  types.foreach { t =>
    jsValues.keySet.diff(t.keySet) match {
      case d if d.nonEmpty =>
        throw new Exception(s"missing type definition(s) ${d.mkString(",")}")
      case _ => ()
    }
  }

  def deserialize(dxFileDescCache: DxFileDescCache = DxFileDescCache.empty,
                  dxApi: DxApi = DxApi.get): Map[String, Value] = {
    val deserializer = ParameterLinkDeserializer(dxFileDescCache, dxApi)
    jsValues.map {
      case (name, jsv) if types.isDefined =>
        val irType = types.get(name)
        name -> deserializer.deserializeInputWithType(jsv, irType, name)
      case (name, jsv) =>
        name -> deserializer.deserializeInput(jsv)
    }
  }

  def toJson: JsObject = {
    val (typesField, allDefinitions, jsDefinitions) = types match {
      case Some(types) if types.nonEmpty =>
        val (jsTypes, allDefinitions, jsDefinitions) =
          types.foldLeft(SortedMap.empty[String, JsValue],
                         definitions.getOrElse(Map.empty),
                         SortedMap.empty[String, JsValue]) {
            case ((typeAccu, defAccu, defJsAccu), (name, t: TSchema)) if defAccu.contains(t.name) =>
              (typeAccu + (name -> JsString(t.name)), defAccu, defJsAccu)
            case ((typeAccu, defAccu, defJsAccu), (name, t: TSchema)) =>
              (typeAccu + (name -> JsString(t.name)), defAccu + (t.name -> t), defJsAccu)
            case ((typeAccu, defAccu, defJsAccu), (name, t)) =>
              val (jsType, newDefJs) = TypeSerde.serialize(t, defJsAccu)
              (typeAccu + (name -> jsType), defAccu, newDefJs)
          }
        (Some(Manifest.TypesKey -> JsObject(jsTypes)), allDefinitions, jsDefinitions)
      case _ =>
        (None, definitions.getOrElse(Map.empty), Map.empty[String, JsValue])
    }
    val definitionsField = if (allDefinitions.nonEmpty) {
      Some(Manifest.DefinitionsKey -> JsObject(allDefinitions.foldLeft(jsDefinitions) {
        case (accu, (name, _)) if accu.contains(name) => accu
        case (accu, (name, schema: TSchema)) =>
          val (schemaJs, newJsDefinitions) = TypeSerde.serializeSchema(schema, accu)
          newJsDefinitions + (name -> schemaJs)
      }))
    } else {
      None
    }
    val fields = Vector(
        definitionsField,
        typesField,
        Some(Manifest.ValuesKey -> JsObject(jsValues)),
        id.map(id => Manifest.IdKey -> JsString(id))
    )
    JsObject(fields.flatten.to(TreeMap))
  }
}
