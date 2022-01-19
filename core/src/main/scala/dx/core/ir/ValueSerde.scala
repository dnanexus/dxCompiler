package dx.core.ir

import dx.api.{DxApi, DxFile, DxFileDescCache}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.CollectionUtils.IterableOnceExtensions
import dx.util.protocols.{DxArchiveFolderSource, DxFileSource, DxFolderSource}
import dx.util.{FileSourceResolver, JsUtils, LocalFileSource}
import spray.json._

import scala.collection.immutable.SeqMap

object ValueSerde extends DefaultJsonProtocol {
  val WrappedValueKey = "wrapped___"

  case class ValueSerdeException(message: String) extends Exception(message)

  /** Serialize a PathValue (VFile or VDirectory). If `fileResolver` is specified, file URIs are
    * resolved and serialized as dx link objects.
    * @param path
    *   PathValue to serialize
    * @param fileResolver
    *   FileSourceResolver
    * @param pathsAsObjects
    *   whether PathValues should be serialized as JsObject or JsString
    * @return
    */
  def serializePath(path: PathValue,
                    fileResolver: Option[FileSourceResolver] = None,
                    pathsAsObjects: Boolean = false
  ): JsValue = {
    def serializeFileUri(uri: String): JsValue = {
      fileResolver
        .map { res =>
          res.resolve(uri) match {
            case dxFile: DxFileSource       => dxFile.dxFile.asJson
            case localFile: LocalFileSource => JsString(localFile.originalPath.toString)
            case other =>
              throw new RuntimeException(s"Unsupported file source ${other}")
          }
        }
        .getOrElse(JsString(uri))
    }
    def serializeFolderUri(uri: String): JsValue = {
      fileResolver
        .map { res =>
          res.resolveDirectory(uri) match {
            case dxFolder: DxFolderSource         => JsString(dxFolder.address)
            case dxArchive: DxArchiveFolderSource => dxArchive.dxFileSource.dxFile.asJson
            case localFile: LocalFileSource       => JsString(localFile.originalPath.toString)
            case other =>
              throw new RuntimeException(s"Unsupported file source ${other}")
          }
        }
        .getOrElse(JsString(uri))
    }
    def inner(innerPath: PathValue): JsValue = {
      innerPath match {
        case VFile(uri, None, None, None, None, Vector(), None) if !pathsAsObjects =>
          serializeFileUri(uri)
        case VFolder(uri, None, None) if !pathsAsObjects =>
          serializeFolderUri(uri)
        case f: VFile if pathsAsObjects =>
          JsObject(
              Vector(
                  Some("type" -> JsString("File")),
                  Some("uri" -> serializeFileUri(f.uri)),
                  f.basename.map("basename" -> JsString(_)),
                  f.contents.map("contents" -> JsString(_)),
                  f.checksum.map("checksum" -> JsString(_)),
                  f.size.map("size" -> JsNumber(_)),
                  Option.when(f.secondaryFiles.nonEmpty)(
                      "secondaryFiles" -> JsArray(f.secondaryFiles.map(inner))
                  ),
                  f.format.map(f => "format" -> JsString(f))
              ).flatten.toMap
          )
        case VFolder(uri, basename, listing) if pathsAsObjects =>
          JsObject(
              Vector(
                  Some("type" -> JsString("Folder")),
                  Some("uri" -> serializeFolderUri(uri)),
                  basename.map("basename" -> JsString(_)),
                  listing.map(l => "listing" -> JsArray(l.map(inner)))
              ).flatten.toMap
          )
        case VListing(basename, listing) =>
          JsObject(
              "type" -> JsString("Listing"),
              "basename" -> JsString(basename),
              "listing" -> JsArray(listing.map(inner))
          )
        case _: PathValue =>
          throw new Exception(s"cannot serialize ${innerPath} with pathsAsObjects=false")
        case _ =>
          throw new Exception(s"not a PathValue ${innerPath}")
      }
    }
    inner(path)
  }

  /** Serializes a Value to JSON.
    * @param value
    *   the Value to serialize
    * @param handler
    *   an optional function to perform special handling of certain values. If Right(jsValue) is
    *   returned, then jsValue is the result of the transformation. If Left(newValue) is returned,
    *   then newValue is transformed according to the default rules.
    * @param pathsAsObjects
    *   whether to serialize path types (File and Directory) as Objects rather than Strings.
    * @return
    */
  def serialize(value: Value,
                handler: Option[Value => Either[Value, JsValue]] = None,
                fileResolver: Option[FileSourceResolver] = None,
                pathsAsObjects: Boolean = false
  ): JsValue = {
    def inner(innerValue: Value): JsValue = {
      val v = handler.map(_(innerValue)) match {
        case Some(Right(result))  => return result
        case Some(Left(newValue)) => newValue
        case None                 => innerValue
      }
      v match {
        case VNull         => JsNull
        case VBoolean(b)   => JsBoolean(b)
        case VInt(i)       => JsNumber(i)
        case VFloat(f)     => JsNumber(f)
        case VString(s)    => JsString(s)
        case p: PathValue  => serializePath(p, fileResolver, pathsAsObjects)
        case VArray(items) => JsArray(items.map(inner))
        case VHash(fields) => JsObject(fields.view.mapValues(inner).toMap)
      }
    }
    inner(value)
  }

  /** Wrap value for a multi-type (including `Any`)
    */
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
      handler: Option[(Value, Type) => Either[Value, JsValue]] = None,
      fileResolver: Option[FileSourceResolver] = None,
      pathsAsObjects: Boolean = false
  ): JsValue = {
    lazy val handlerWrapper = handler.map(h => (v: Value) => h(v, TMulti.Any))
    def inner(innerValue: Value, innerType: Type): JsValue = {
      val v = handler.map(_(innerValue, innerType)) match {
        case Some(Right(result))  => return result
        case Some(Left(newValue)) => newValue
        case None                 => innerValue
      }
      (innerType, v) match {
        case (_: TOptional, VNull)             => JsNull
        case (TOptional(t), _)                 => inner(v, t)
        case (t: TMulti, _)                    => innerMulti(v, t)
        case (TBoolean, VBoolean(b))           => JsBoolean(b)
        case (TInt, VInt(i))                   => JsNumber(i)
        case (TInt, VFloat(i)) if i.isValidInt => JsNumber(i.intValue())
        case (TFloat, VFloat(f))               => JsNumber(f)
        case (TFloat, VInt(i))                 => JsNumber(i.floatValue())
        case (TString, VString(s))             => JsString(s)
        case (TString, f: VFile)               => JsString(f.uri)
        case (TString, f: VFolder)             => JsString(f.uri)
        case (TFile, f: VFile)                 => serializePath(f, fileResolver, pathsAsObjects)
        case (TDirectory, d: DirectoryValue)   => serializePath(d, fileResolver, pathsAsObjects)
        case (TArray(_, true), VArray(items)) if items.isEmpty =>
          throw ValueSerdeException(s"empty value for non-empty array type ${innerType}")
        case (TArray(itemType, _), VArray(items)) =>
          JsArray(items.map(inner(_, itemType)))
        case (TSchema(schemaName, fieldTypes), VHash(fields)) =>
          if (!fields.keySet.subsetOf(fieldTypes.keySet)) {
            throw ValueSerdeException(
                s"${fields} contains invalid field(s) (schema ${schemaName})"
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
          .collectFirstDefined { t =>
            try {
              Some(inner(innerValue, t))
            } catch {
              case _: Throwable => None
            }
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
    values.map { case (name, value) =>
      name -> serialize(value)
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

  def isPathObject(jsv: JsValue): Boolean = {
    jsv match {
      case JsObject(fields) =>
        fields.get("type") match {
          case Some(JsString("File" | "Folder" | "Archive")) if fields.contains("uri") => true
          case Some(JsString("Listing")) if fields.contains("listing")                 => true
          case _                                                                       => false
        }
      case _ => false
    }
  }

  def deserializeDxFileUri(uri: JsValue,
                           dxApi: DxApi = DxApi.get,
                           dxFileDescCache: Option[DxFileDescCache] = None
  ): String = {
    uri match {
      case JsString(s) => s
      case obj: JsObject if DxFile.isLinkJson(obj) =>
        val dxFile = DxFile.fromJson(dxApi, obj)
        dxFileDescCache.map(_.updateFileFromCache(dxFile)).getOrElse(dxFile).asUri
      case _ =>
        throw new Exception(s"invalid URI value ${uri}")
    }
  }

  def deserializePathObject(jsv: JsValue,
                            dxApi: DxApi = DxApi.get,
                            dxFileDescCache: Option[DxFileDescCache] = None
  ): PathValue = {
    def deserializeListing(fields: Map[String, JsValue]): Option[Vector[PathValue]] = {
      JsUtils
        .getOptionalValues(fields, "listing")
        .map {
          _.map {
            case obj: JsObject => deserializePathObject(obj)
            case other         => throw new Exception(s"invalid path ${other}")
          }
        }
    }
    def inner(innerValue: JsValue): PathValue = {
      innerValue match {
        case JsObject(fields) =>
          fields.get("type") match {
            case Some(JsString("File")) =>
              VFile(
                  deserializeDxFileUri(fields("uri"), dxApi, dxFileDescCache),
                  JsUtils.getOptionalString(fields, "basename"),
                  JsUtils.getOptionalString(fields, "contents"),
                  JsUtils.getOptionalString(fields, "checksum"),
                  JsUtils.getOptionalLong(fields, "size"),
                  JsUtils
                    .getOptionalValues(fields, "secondaryFiles")
                    .map {
                      _.map {
                        case obj: JsObject => deserializePathObject(obj)
                        case other         => throw new Exception(s"invalid path ${other}")
                      }
                    }
                    .getOrElse(Vector.empty),
                  JsUtils.getOptionalString(fields, "format")
              )
            case Some(JsString("Folder")) =>
              VFolder(
                  deserializeDxFileUri(fields("uri"), dxApi, dxFileDescCache),
                  JsUtils.getOptionalString(fields, "basename"),
                  deserializeListing(fields)
              )
            case Some(JsString("Listing")) =>
              VListing(
                  JsUtils.getString(fields, "basename"),
                  deserializeListing(fields).getOrElse(Vector.empty)
              )
            case _ =>
              throw new Exception(s"invalid serialized PathValue ${fields}")
          }
        case _ => throw new Exception("")
      }
    }
    inner(jsv)
  }

  /** Deserializes a JsValue to a Value, in the absence of type information.
    *
    * @param jsValue
    *   the JsValue
    * @param handler
    *   an optional function to perform special handling of certain values. If Right(value) is
    *   returned, then value is the result of the transformation. If Left(newJsValue) is returned,
    *   then newJsValue is transformed according to the default rules.
    * @return
    */
  def deserialize(jsValue: JsValue,
                  handler: Option[JsValue => Either[JsValue, Value]] = None,
                  dxApi: DxApi = DxApi.get,
                  dxFileDescCache: Option[DxFileDescCache] = None
  ): Value = {
    def inner(innerValue: JsValue): Value = {
      val v = handler.map(_(innerValue)) match {
        case Some(Right(result))    => return result
        case Some(Left(newJsValue)) => newJsValue
        case None                   => innerValue
      }
      v match {
        case _ if isWrappedValue(v)                         => inner(unwrapValue(v))
        case JsNull                                         => VNull
        case JsTrue                                         => VBoolean(true)
        case JsFalse                                        => VBoolean(false)
        case JsNumber(value) if value.isValidLong           => VInt(value.toLongExact)
        case JsNumber(value)                                => VFloat(value.toDouble)
        case JsString(s) if DxFileSource.isDxFileUri(s)     => VFile(s)
        case JsString(s) if DxFolderSource.isDxFolderUri(s) => VFolder(s)
        case JsString(s)                                    => VString(s)
        case JsArray(items)                                 => VArray(items.map(x => inner(x)))
        case obj: JsObject if isPathObject(obj) =>
          deserializePathObject(obj, dxApi, dxFileDescCache)
        case obj: JsObject if DxFile.isLinkJson(obj) =>
          VFile(deserializeDxFileUri(obj, dxApi, dxFileDescCache))
        case JsObject(fields) => VHash(fields.view.mapValues(inner).to(SeqMap))
      }
    }
    inner(jsValue)
  }

  /** Deserializes a JsValue to a Value of the specified type.
    * @param jsValue
    *   the JsValue
    * @param t
    *   the Type
    * @param handler
    *   an optional function to perform special handling of certain values. If Right(value) is
    *   returned, then value is the result of the transformation. If Left(newJsValue) is returned,
    *   then newJsValue is transformed according to the default rules.
    * @return
    */
  def deserializeWithType(
      jsValue: JsValue,
      t: Type,
      name: String = "",
      handler: Option[(JsValue, Type) => Either[JsValue, Value]] = None,
      dxApi: DxApi = DxApi.get,
      dxFileDescCache: Option[DxFileDescCache] = None
  ): Value = {
    def inner(innerValue: JsValue, innerType: Type, innerName: String): Value = {
      val v = handler.map(_(innerValue, innerType)) match {
        case Some(Right(result))    => return result
        case Some(Left(newJsValue)) => newJsValue
        case None                   => innerValue
      }
      (innerType, v) match {
        case (TOptional(_), JsNull) => VNull
        case (TOptional(t), _)      => inner(v, t, innerName)
        case (any: TMulti, _) if any.bounds.isEmpty && isWrappedValue(v) =>
          inner(unwrapValue(v), any, innerName)
        case (TMulti(bounds), _) if bounds.isEmpty => deserialize(v)
        case (TMulti(bounds), _) =>
          val unwrappedValue = unwrapValue(v)
          bounds.iterator
            .collectFirstDefined { t =>
              try {
                Some(inner(unwrappedValue, t, innerName))
              } catch {
                case _: ValueSerdeException => None
              }
            }
            .getOrElse(
                throw ValueSerdeException(
                    s"${innerName} value ${unwrappedValue} does not match any of ${bounds}"
                )
            )
        case (TBoolean, JsBoolean(b))                     => VBoolean(b.booleanValue)
        case (TInt, JsNumber(value)) if value.isValidLong => VInt(value.toLongExact)
        case (TFloat, JsNumber(value))                    => VFloat(value.toDouble)
        case (TString, JsString(s))                       => VString(s)
        case (TFile, JsString(uri))                       => VFile(uri)
        case (TDirectory, JsString(uri)) if DxFolderSource.isDxFolderUri(uri) => VFolder(uri)
        case (TFile | TDirectory, obj: JsObject) if isPathObject(obj) =>
          deserializePathObject(obj, dxApi, dxFileDescCache)
        case (TFile, obj: JsObject) if DxFile.isLinkJson(obj) =>
          VFile(deserializeDxFileUri(obj, dxApi, dxFileDescCache))
        case (TArray(_, true), JsArray(items)) if items.isEmpty =>
          throw ValueSerdeException(
              s"Cannot convert ${innerName} empty array to non-empty type ${innerType}"
          )
        case (TArray(t, _), JsArray(items)) =>
          VArray(items.zipWithIndex.map { case (x, index) =>
            inner(x, t, s"${innerName}[${index}]")
          })
        case (TSchema(name, fieldTypes), JsObject(fields)) =>
          // ensure 1) fields keys are a subset of typeTypes keys, 2) fields
          // values are convertible to the corresponding types, and 3) any keys
          // in fieldTypes that do not appear in fields are optional
          val keys1 = fields.keySet
          val keys2 = fieldTypes.keySet
          if (!keys1.subsetOf(keys2)) {
            throw ValueSerdeException(
                s"keys (${keys1}) have members that do not appear in struct ${name}"
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
            case (name, t) if fields.contains(name) =>
              name -> inner(fields(name), t, s"${innerName}.${name}")
          })
        case (THash, JsObject(fields)) =>
          VHash(
              fields
                .map { case (key, value) =>
                  key -> deserialize(value)
                }
                .to(SeqMap)
          )
        case (TEnum(symbols), JsString(s)) if symbols.contains(s) =>
          VString(s)
        case _ =>
          throw ValueSerdeException(
              s"cannot deserialize ${innerName} value ${innerValue} as type ${innerType}"
          )
      }
    }
    inner(jsValue, t, name)
  }

  def deserializeMap(m: Map[String, JsValue]): Map[String, Value] = {
    m.map { case (k, v) =>
      k -> deserialize(v)
    }
  }

  // support automatic conversion to/from JsValue
  implicit val valueFormat: RootJsonFormat[Value] = new RootJsonFormat[Value] {
    override def read(jsv: JsValue): Value = deserialize(jsv)
    override def write(value: Value): JsValue = serialize(value)
  }

  def toString(value: Value, verbose: Boolean = false, indent: Int = 0): String = {
    val prefix = " " * 2 * indent
    def pathToString(p: PathValue): String = {
      p match {
        case f: VFile if verbose =>
          Vector(
              Some(s"${prefix}${f.uri}"),
              f.basename.map(bn => s"basename: ${bn}"),
              f.contents.map(c => s"contents: ${c}"),
              Option.when(f.secondaryFiles.nonEmpty)(
                  s"secondaryFiles:\n${f.secondaryFiles.map(toString(_, verbose = true, indent = indent + 2))}"
              )
          ).flatten.mkString(s"\n${prefix}  ")
        case f: VFile => s"${prefix}${f.uri}"
        case f: VFolder if verbose =>
          Vector(
              Some(s"${prefix}${f.uri}"),
              f.basename.map(bn => s"basename: ${bn}"),
              Option.when(f.listing.nonEmpty)(
                  s"listing:\n${f.listing.get.map(toString(_, verbose = true, indent = indent + 2))}"
              )
          ).flatten.mkString(s"\n${prefix}  ")
        case f: VFolder => s"${prefix}${f.uri}"
        case l: VListing if verbose =>
          s"${prefix}${l.basename}\n${prefix}  listing:\n${l.items
            .map(toString(_, verbose = true, indent = indent + 2))}"
        case l: VListing => s"${prefix}${l.basename}"
      }
    }
    value match {
      case p: PathValue => pathToString(p)
      case VInt(i)      => s"${prefix}${i.toString}"
      case VFloat(f)    => s"${prefix}${f.toString}"
      case VString(s)   => s"${prefix}${s}"
      case VBoolean(b)  => s"${prefix}${b.toString}"
      case VNull        => s"${prefix}null"
      case VArray(items) if verbose =>
        s"${prefix}[\n${items.map(toString(_, verbose = true, indent = indent + 1)).mkString("\n")}\n${prefix}]"
      case VArray(items) => s"[${items.map(toString(_)).mkString(",")}]"
      case VHash(fields) if verbose =>
        s"${prefix}{\n${fields
          .map { case (k, v) => s"${k}: ${toString(v, verbose = true, indent = indent + 1)}" }
          .mkString("\n")}\n${prefix}}"
      case VHash(fields) => s"{${fields.map { case (k, v) => s"${k}: ${toString(v)}" }}}"
      case _             => throw new Exception(s"unexpected value ${value}")
    }
  }
}
