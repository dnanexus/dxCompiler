package dx.core.ir

import dx.api.{DxApi, DxFileDescCache}
import dx.core.ir.Type._
import dx.util.JsUtils
import spray.json._

import java.nio.file.Path
import scala.collection.immutable.{SortedMap, TreeMap}

object Manifest {
  val IdKey = "id"
  val EncodedKey = "encoded"
  val ValuesKey = "values"
  val TypesKey = "types"
  val DefinitionsKey = "definitions"

  val AllKeys = Set(IdKey, EncodedKey, ValuesKey, TypesKey, DefinitionsKey)

  def keysToNames[T](values: Map[String, T],
                     dxNameFactory: DxNameFactory,
                     encoded: Boolean = false
  ): Map[DxName, T] = {
    values.map {
      case (name, value) if encoded => dxNameFactory.fromEncodedName(name) -> value
      case (name, value)            => dxNameFactory.fromDecodedName(name) -> value
    }
  }

  def parseFullManifest(jsValue: JsValue,
                        dxNameFactory: DxNameFactory,
                        irTypes: Option[Map[DxName, Type]] = None,
                        typeAliases: Map[String, Type] = Map.empty
  ): Manifest = {
    val fields = jsValue.asJsObject.fields
    val id = fields.get(Manifest.IdKey) match {
      case Some(JsString(id)) => Some(id)
      case None               => None
      case other              => throw new Exception(s"invalid manifest ID ${other}")
    }
    val encoded = fields.get(Manifest.EncodedKey) match {
      case Some(JsBoolean(encoded)) => encoded
      case None                     => false
      case other                    => throw new Exception(s"invalid encoded value ${other}")
    }
    val jsValues: Map[DxName, JsValue] = fields.get(Manifest.ValuesKey) match {
      case Some(JsObject(values)) => keysToNames(values, dxNameFactory, encoded)
      case other                  => throw new Exception(s"invalid manifest values ${other}")
    }
    val jsDefinitions: Map[String, JsValue] = fields.get(Manifest.DefinitionsKey) match {
      case Some(JsObject(fields)) => fields
      case None                   => Map.empty[String, JsValue]
      case other                  => throw new Exception(s"invalid 'definitions' ${other}")
    }
    val (types, definitions) = fields.get(Manifest.TypesKey) match {
      case Some(JsObject(typeDefs)) =>
        val (types, definitions) =
          TypeSerde.deserializeMap(typeDefs, typeAliases, jsDefinitions) match {
            case (types, definitions) => (keysToNames(types, dxNameFactory, encoded), definitions)
          }
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
              case _ => ()
            }
        }
        (Some(types), Some(definitions))
      case None => (irTypes, None)
      case _ =>
        throw new Exception(s"invalid manifest ${jsValue}")
    }
    Manifest(jsValues, types, definitions, id)
  }

  def parse(jsValue: JsValue,
            dxNameFactory: DxNameFactory,
            types: Option[Map[DxName, Type]] = None,
            typeAliases: Map[String, Type] = Map.empty
  ): Manifest = {
    jsValue match {
      case JsObject(fields)
          if fields.contains(Manifest.ValuesKey) &&
            fields.keySet.diff(Manifest.AllKeys).isEmpty =>
        parseFullManifest(jsValue, dxNameFactory, types, typeAliases)
      case _ =>
        Manifest(
            keysToNames(jsValue.asJsObject.fields, dxNameFactory),
            types,
            types.map(t => Type.collectSchemas(t.values.toVector))
        )
    }
  }

  def parseFile(path: Path,
                dxNameFactory: DxNameFactory,
                types: Option[Map[DxName, Type]] = None,
                typeAliases: Map[String, Type] = Map.empty
  ): Manifest = {
    parse(JsUtils.jsFromFile(path), dxNameFactory, types, typeAliases)
  }
}

case class Manifest(jsValues: Map[DxName, JsValue],
                    types: Option[Map[DxName, Type]] = None,
                    definitions: Option[Map[String, Type]] = None,
                    id: Option[String] = None
) {
  types.foreach { t =>
    jsValues.keySet.diff(t.keySet) match {
      case d if d.nonEmpty =>
        throw new Exception(s"missing type definition(s) ${d.mkString(",")}")
      case _ => ()
    }
  }

  def deserialize(dxFileDescCache: DxFileDescCache = DxFileDescCache.empty,
                  dxApi: DxApi = DxApi.get
  ): Map[DxName, Value] = {
    val deserializer = ParameterLinkDeserializer(dxFileDescCache, dxApi)
    jsValues.map {
      case (dxName, jsv) if types.isDefined =>
        val irType = types.get(dxName)
        dxName -> deserializer.deserializeInputWithType(jsv, irType, dxName.decoded)
      case (dxName, jsv) =>
        dxName -> deserializer.deserializeInput(jsv)
    }
  }

  def toJson(encodeNames: Boolean = false): JsObject = {
    def dxNameToString(dxName: DxName): String = {
      if (encodeNames) {
        dxName.encoded
      } else {
        dxName.decoded
      }
    }
    val (typesField, allDefinitions, jsDefinitions) = types match {
      case Some(types) if types.nonEmpty =>
        val (jsTypes, allDefinitions, jsDefinitions) =
          types.foldLeft(SortedMap.empty[String, JsValue],
                         definitions.getOrElse(Map.empty),
                         SortedMap.empty[String, JsValue]
          ) {
            case ((typeAccu, defAccu, defJsAccu), (dxName, t: TSchema))
                if defAccu.contains(t.name) =>
              (typeAccu + (dxNameToString(dxName) -> JsString(t.name)), defAccu, defJsAccu)
            case ((typeAccu, defAccu, defJsAccu), (dxName, t: TSchema)) =>
              (typeAccu + (dxNameToString(dxName) -> JsString(t.name)),
               defAccu + (t.name -> t),
               defJsAccu
              )
            case ((typeAccu, defAccu, defJsAccu), (dxName, t)) =>
              val (jsType, newDefJs) = TypeSerde.serialize(t, defJsAccu)
              (typeAccu + (dxNameToString(dxName) -> jsType), defAccu, newDefJs)
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
        case (_, (_, other)) => throw new Exception(s"invalid schema ${other}")
      }))
    } else {
      None
    }
    val fields = Vector(
        definitionsField,
        typesField,
        Some(Manifest.ValuesKey -> JsObject(jsValues.map { case (dxName, jsv) =>
          dxNameToString(dxName) -> jsv
        })),
        id.map(id => Manifest.IdKey -> JsString(id)),
        Some(Manifest.EncodedKey -> JsBoolean(encodeNames))
    )
    JsObject(fields.flatten.to(TreeMap))
  }
}
