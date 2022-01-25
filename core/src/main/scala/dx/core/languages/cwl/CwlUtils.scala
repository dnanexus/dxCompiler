package dx.core.languages.cwl

import dx.api.DxPath
import dx.core.io.DxWorkerPaths
import dx.core.ir.{DxName, Type, Value}
import dx.core.ir.Type._
import dx.core.ir.Value.{
  VArray,
  VBoolean,
  VFile,
  VFloat,
  VFolder,
  VHash,
  VInt,
  VListing,
  VNull,
  VString,
  DirectoryValue => IRDirectoryValue,
  PathValue => IRPathValue
}
import dx.cwl._
import dx.util.CollectionUtils.IterableOnceExtensions
import dx.util.{FileSourceResolver, LocalFileSource}
import spray.json._

import java.nio.file.{Path, Paths}
import java.util.UUID
import scala.annotation.tailrec
import scala.collection.immutable.SeqMap

object CwlUtils {

  /**
    * CWL has a "null" type. When used alone, it means the field
    * only accepts a null value. This is not able to be represented
    * natively, so instead we use an optional schema type.
    */
  val NullSchema: Type = TOptional(TSchema("null", SeqMap.empty))

  def isNullSchema(t: TSchema): Boolean = {
    t.name == "null"
  }

  def isNullValue(fields: Map[String, _]): Boolean = {
    fields.isEmpty
  }

  def toIRSchema(cwlRecord: CwlRecord): TSchema = {
    if (!cwlRecord.hasName) {
      throw new Exception(s"cannot convert schema without name ${cwlRecord}")
    }
    TSchema(
        cwlRecord.name,
        cwlRecord.fields.map {
          case (key, value) => key -> toIRType(value.cwlType)
        }
    )
  }

  def toIRType(cwlType: CwlType): Type = {
    cwlType match {
      case CwlOptional(t)  => Type.ensureOptional(toIRType(t))
      case CwlMulti(types) => TMulti(types.map(toIRType))
      case CwlAny          => TMulti.Any
      case CwlNull         => NullSchema
      case CwlBoolean      => TBoolean
      case CwlInt          => TInt
      case CwlLong         => TInt
      case CwlDouble       => TFloat
      case CwlFloat        => TFloat
      case CwlString       => TString
      case CwlFile         => TFile
      case CwlDirectory    => TDirectory
      case a: CwlArray     => TArray(toIRType(a.itemType))
      case e: CwlEnum      => TEnum(e.symbolNames)
      case r: CwlRecord if r.hasName =>
        toIRSchema(r)
      case _ =>
        throw new Exception(s"Cannot convert CWL type ${cwlType} to IR")
    }
  }

  def toIRPath(path: PathValue): IRPathValue = {
    path match {
      case f @ FileValue(location,
                         path,
                         basename,
                         _,
                         _,
                         _,
                         checksum,
                         size,
                         secondaryFiles,
                         format,
                         contents) =>
        VFile(
            location
              .orElse(path)
              .orElse(
                  Option.when(contents.isDefined)(basename.getOrElse(UUID.randomUUID().toString))
              )
              .getOrElse(
                  throw new Exception(s"FileValue does not have 'location' or 'contents' ${f}")
              ),
            basename,
            contents,
            checksum,
            size,
            secondaryFiles.map(toIRPath),
            format
        )
      case d @ DirectoryValue(location, path, basename, listing) =>
        location.orElse(path) match {
          case Some(uri) =>
            VFolder(uri, basename, Option.unless(listing.isEmpty)(listing.map(toIRPath)))
          case None if listing.nonEmpty =>
            VListing(basename.get, listing.map(toIRPath))
          case _ =>
            throw new Exception(
                s"DirectoryValue does not have 'location', 'path', or 'listing' ${d}"
            )
        }
    }
  }

  private def toIRValue(cwlValue: CwlValue): (Type, Value) = {
    cwlValue match {
      case BooleanValue(b)   => (TBoolean, VBoolean(value = b))
      case IntValue(i)       => (TInt, VInt(i))
      case LongValue(l)      => (TInt, VInt(l))
      case FloatValue(f)     => (TFloat, VFloat(f))
      case DoubleValue(d)    => (TFloat, VFloat(d))
      case StringValue(s)    => (TString, VString(s))
      case f: FileValue      => (TFile, toIRPath(f))
      case d: DirectoryValue => (TDirectory, toIRPath(d))
      case NullValue         =>
        // it is possible that the type is null (i.e. NullSchema) but much more likely CwlOptional
        (TMulti.Any, VNull)
      case ArrayValue(items) =>
        val (itemTypes, itemValues, optional) =
          items.foldLeft(Set.empty[Type], Vector.empty[Value], false) {
            case ((types, values, optional), item) =>
              item match {
                case NullValue => (types, values :+ VNull, true)
                case _ =>
                  val (t, v) = toIRValue(item)
                  (types + t, values :+ v, optional)
              }
          }
        val itemType = itemTypes.toVector match {
          case Vector()              => TMulti.Any
          case Vector(t) if optional => ensureOptional(t)
          case Vector(t)             => t
          case bounds if optional    => Type.merge(bounds.map(ensureOptional(_)))
          case bounds                => Type.merge(bounds)
        }
        (TArray(itemType), VArray(itemValues))
      case ObjectValue(m) =>
        (THash,
         VHash(
             m.map {
                 case (key, NullValue) => key -> VNull
                 case (key, value) =>
                   val (_, v) = toIRValue(value)
                   key -> v
               }
               .to(SeqMap)
         ))
      case _ => throw new Exception(s"Invalid CWL value ${cwlValue})")
    }
  }

  def toIRValue(cwlValue: CwlValue, cwlType: CwlType): (Type, Value) = {
    (cwlType, cwlValue) match {
      case (CwlAny, _) =>
        // TODO: do we need to keep the type as Any or can we return the inferred type?
        (TMulti.Any, toIRValue(cwlValue)._2)
      case (CwlNull, NullValue) => (NullSchema, VNull)
      case (CwlNull, _) =>
        throw new Exception("type 'null' only accepts a value of 'null'")
      case (CwlOptional(_), NullValue) => (toIRType(cwlType), VNull)
      case (CwlOptional(t), _) =>
        val (irType, irValue) = toIRValue(cwlValue, t)
        (Type.ensureOptional(irType), irValue)
      case (CwlMulti(types), _) =>
        types
          .collectFirstDefined {
            case CwlAny => None
            case t =>
              try {
                Some(toIRValue(cwlValue, t))
              } catch {
                case _: Throwable => None
              }
          }
          .getOrElse(
              if (types.contains(CwlAny)) {
                toIRValue(cwlValue)
              } else {
                throw new Exception(
                    s"cannot translate ${cwlValue} as any of ${types.mkString("\n")}"
                )
              }
          )
      case (CwlBoolean, BooleanValue(b))     => (TBoolean, VBoolean(b))
      case (CwlInt, IntValue(i))             => (TInt, VInt(i))
      case (CwlLong, LongValue(l))           => (TInt, VInt(l))
      case (CwlFloat, FloatValue(f))         => (TFloat, VFloat(f))
      case (CwlDouble, DoubleValue(d))       => (TFloat, VFloat(d))
      case (t: CwlNumber, n: NumericValue)   => toIRValue(n.coerceTo(t)._2, t)
      case (CwlString, StringValue(s))       => (TString, VString(s))
      case (CwlFile, f: FileValue)           => (TFile, toIRPath(f))
      case (CwlFile, StringValue(s))         => (TFile, VFile(s))
      case (CwlDirectory, d: DirectoryValue) => (TDirectory, toIRPath(d))
      case (CwlDirectory, StringValue(s))    => (TDirectory, VFolder(s))
      case (array: CwlArray, ArrayValue(items)) =>
        val (irItems, optional) =
          items.foldLeft(Vector.empty[Value], false) {
            case ((irItems, _), NullValue) => (irItems :+ VNull, true)
            case ((irItems, optional), i) =>
              val (t, v) = toIRValue(i, array.itemType)
              (irItems :+ v, optional || Type.isOptional(t))
          }
        val irItemType = toIRType(array.itemType)
        val irType = TArray(if (optional) {
          Type.ensureOptional(irItemType)
        } else {
          irItemType
        })
        (irType, VArray(irItems))
      case (record: CwlRecord, ObjectValue(fields)) =>
        // generate a random name for anonymous schemas
        val name = if (record.hasName) record.name else UUID.randomUUID().toString
        val (types, values) =
          fields.foldLeft(SeqMap.empty[String, Type], SeqMap.empty[String, Value]) {
            case ((types, values), (name, value)) if record.fields.contains(name) =>
              val (irType, irValue) = toIRValue(value, record.fields(name).cwlType)
              (types + (name -> irType), values + (name -> irValue))
            case (name, _) =>
              throw new Exception(s"invalid field ${name}")
          }
        (TSchema(name, types), VHash(values))
      case (enum: CwlEnum, StringValue(s)) if enum.symbolNames.contains(s) =>
        (TEnum(enum.symbolNames), VString(s))
      case _ =>
        throw new Exception(
            s"""Cannot convert CWL value ${prettyFormatValue(cwlValue, withType = true)}
               |to type ${prettyFormatType(cwlType)})""".stripMargin.replaceAll("\n", " ")
        )
    }
  }

  def toIR(cwl: Map[DxName, (CwlType, CwlValue)]): Map[DxName, (Type, Value)] = {
    cwl.map {
      case (name, (cwlType, cwlValue)) => name -> toIRValue(cwlValue, cwlType)
    }
  }

  def fromIRType(irType: Type,
                 typeAliases: Map[String, CwlSchema] = Map.empty,
                 isInput: Boolean): CwlType = {
    def inner(innerType: Type): CwlType = {
      innerType match {
        case TOptional(t)  => CwlOptional(inner(t))
        case TMulti.Any    => CwlAny
        case TMulti(types) => CwlMulti(types.map(inner))
        case NullSchema    => CwlNull
        case TBoolean      => CwlBoolean
        case TInt          => CwlLong
        case TFloat        => CwlDouble
        case TString       => CwlString
        case TFile         => CwlFile
        case TDirectory    => CwlDirectory
        case TArray(t, _)  => CwlArray(inner(t))
        case TSchema(name, _) if typeAliases.contains(name) =>
          typeAliases(name)
        case TSchema(name, fields) if isInput =>
          CwlInputRecord(fields.map {
            case (name, t) => name -> CwlInputRecordField(name, inner(t))
          }, Some(Identifier(namespace = None, frag = name)))
        case TSchema(name, members) =>
          CwlOutputRecord(members.map {
            case (name, t) => name -> CwlOutputRecordField(name, inner(t))
          }, Some(Identifier(namespace = None, frag = name)))
        case TEnum(symbols) => CwlEnum(symbols)
        case _ =>
          throw new Exception(s"Cannot convert IR type ${irType} to CWL")
      }
    }
    inner(irType)
  }

  def fromIRPath(path: IRPathValue): PathValue = {
    path match {
      case f: VFile =>
        FileValue(Some(f.uri),
                  basename = f.basename,
                  checksum = f.checksum,
                  size = f.size,
                  contents = f.contents,
                  secondaryFiles = f.secondaryFiles.map(fromIRPath),
                  format = f.format)
      case VFolder(uri, basename, listing) =>
        DirectoryValue(location = Some(uri),
                       basename = basename,
                       listing = listing.map(_.map(fromIRPath)).getOrElse(Vector.empty[PathValue]))
      case VListing(basename, listing) =>
        DirectoryValue(basename = Some(basename), listing = listing.map(fromIRPath))
    }
  }

  def fromIRValue(
      value: Value,
      name: Option[String],
      isInput: Boolean,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      resolveLocalPaths: Boolean = false
  ): (CwlType, CwlValue) = {
    def inner(innerValue: Value, innerName: Option[String]): (CwlType, CwlValue) = {
      innerValue match {
        case VNull       => (CwlOptional(CwlAny), NullValue)
        case VBoolean(b) => (CwlBoolean, BooleanValue(b))
        case VInt(i)     => (CwlLong, LongValue(i))
        case VFloat(f)   => (CwlDouble, DoubleValue(f))
        case VString(s) =>
          fileResolver.resolve(s) match {
            case local: LocalFileSource if resolveLocalPaths && local.exists && local.isDirectory =>
              (CwlDirectory, DirectoryValue(s))
            case local: LocalFileSource if resolveLocalPaths && local.exists =>
              (CwlFile, FileValue(s))
            case _: LocalFileSource   => (CwlString, StringValue(s))
            case fs if fs.isDirectory => (CwlDirectory, DirectoryValue(s))
            case _                    => (CwlFile, FileValue(s))
          }
        case f: VFile            => (CwlFile, fromIRPath(f))
        case d: IRDirectoryValue => (CwlDirectory, fromIRPath(d))
        case VArray(array) =>
          val (types, values) = array.zipWithIndex.map {
            case (v, i) => inner(v, innerName.map(n => s"${n}[${i}]"))
          }.unzip
          (CwlArray(CwlType.flatten(types.distinct)), ArrayValue(values))
        case VHash(fields) =>
          val (types, values) = fields.map {
            case (key, value) =>
              val (cwlType, cwlValue) = inner(value, innerName.map(n => s"${n}[${key}]"))
              (key -> cwlType, key -> cwlValue)
          }.unzip
          // create an anonymous record schema
          val schemaType = if (isInput) {
            CwlInputRecord(
                types
                  .map {
                    case (name, t) => name -> CwlInputRecordField(name, t)
                  }
                  .to(SeqMap)
            )
          } else {
            CwlOutputRecord(
                types
                  .map {
                    case (name, t) => name -> CwlOutputRecordField(name, t)
                  }
                  .to(SeqMap)
            )
          }
          (schemaType, ObjectValue(values.to(SeqMap)))
        case _ =>
          throw new Exception(
              s"cannot convert ${name.getOrElse("IR")} value ${value} to WDL value"
          )
      }
    }
    inner(value, name)
  }

  def fromIRValues(
      values: Map[String, Value],
      isInput: Boolean,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      resolveLocalPaths: Boolean = false
  ): Map[String, (CwlType, CwlValue)] = {
    values.map {
      case (name, value) =>
        name -> fromIRValue(value,
                            Some(name),
                            isInput,
                            fileResolver = fileResolver,
                            resolveLocalPaths = resolveLocalPaths)
    }
  }

  def fromIRValueWithType(
      value: Value,
      cwlType: CwlType,
      name: String,
      isInput: Boolean,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      resolveLocalPaths: Boolean = false
  ): (CwlType, CwlValue) = {
    def isPath(path: String): Boolean = {
      fileResolver.resolve(path) match {
        case local: LocalFileSource => resolveLocalPaths && local.exists
        case _                      => true
      }
    }
    @tailrec
    def inner(innerValue: Value, innerType: CwlType, innerName: String): CwlValue = {
      (innerType, innerValue) match {
        case (CwlOptional(_) | CwlNull, VNull) => NullValue
        case (CwlOptional(t), _)               => inner(innerValue, t, innerName)
        case (CwlNull, _) =>
          throw new Exception("type 'null' only accepts a value of 'null'")
        case (CwlBoolean, VBoolean(b))                     => BooleanValue(b)
        case (CwlInt, VInt(i)) if i.isValidInt             => IntValue(i)
        case (CwlLong, VInt(l))                            => LongValue(l)
        case (CwlFloat, VFloat(f))                         => FloatValue(f.toFloat)
        case (CwlFloat, VInt(i))                           => FloatValue(i.toFloat)
        case (CwlDouble, VFloat(f))                        => FloatValue(f)
        case (CwlDouble, VInt(i))                          => FloatValue(i.toDouble)
        case (CwlString, VString(s))                       => StringValue(s)
        case (CwlFile, VString(path)) if isPath(path)      => FileValue(path)
        case (CwlFile, f: VFile)                           => fromIRPath(f)
        case (CwlDirectory, VString(path)) if isPath(path) => DirectoryValue(path)
        case (CwlDirectory, d: IRDirectoryValue)           => fromIRPath(d)
        case (array: CwlArray, VArray(items)) =>
          ArrayValue(items.zipWithIndex.map {
            case (item, i) =>
              innerMulti(item, array.itemType, s"${innerName}[${i}]")._2
          })
        case (record: CwlRecord, VHash(fields)) =>
          // ensure 1) members keys are a subset of memberTypes keys, 2) members
          // values are convertable to the corresponding types, and 3) any keys
          // in memberTypes that do not appear in members are optional
          val keys1 = fields.keySet
          val keys2 = record.fields.keySet
          if (!keys1.subsetOf(keys2)) {
            throw new Exception(
                s"keys (${keys1}) have members that do not appear in struct ${record.name}"
            )
          }
          val missingNonOptional =
            keys1.diff(keys2).map(key => key -> record.fields(key)).filterNot {
              case (_, field) if CwlOptional.isOptional(field.cwlType) => false
              case _                                                   => true
            }
          if (missingNonOptional.nonEmpty) {
            throw new Exception(
                s"struct ${record.name} value is missing non-optional members ${missingNonOptional}"
            )
          }
          ObjectValue(
              fields.map {
                case (key, value) =>
                  key -> innerMulti(value, record.fields(key).cwlType, s"${innerName}[${key}]")._2
              }
          )
        case (enum: CwlEnum, VString(s)) if enum.symbolNames.contains(s) => StringValue(s)
        case _ =>
          throw new Exception(s"cannot translate ${innerValue} to CwlValue of type ${innerType}")
      }
    }
    def innerMulti(innerValue: Value,
                   innerType: CwlType,
                   innerName: String): (CwlType, CwlValue) = {
      innerType match {
        case CwlMulti(types) =>
          types
            .collectFirstDefined {
              case CwlAny => None
              case t =>
                try {
                  Some(t, inner(innerValue, t, innerName))
                } catch {
                  case _: Throwable => None
                }
            }
            .getOrElse(
                if (types.contains(CwlAny)) {
                  fromIRValue(innerValue,
                              Some(innerName),
                              isInput,
                              fileResolver = fileResolver,
                              resolveLocalPaths = resolveLocalPaths)
                } else {
                  throw new Exception(
                      s"""cannot convert ${innerName} ${innerValue} to CWL value of any type 
                         |${types.mkString(",")}""".stripMargin.replaceAll("\n", " ")
                  )
                }
            )
        case CwlOptional(t: CwlMulti) =>
          val (cwlType, cwlValue) = innerMulti(innerValue, t, innerName)
          (CwlOptional.ensureOptional(cwlType), cwlValue)
        case CwlAny =>
          fromIRValue(innerValue,
                      Some(innerName),
                      isInput,
                      fileResolver = fileResolver,
                      resolveLocalPaths = resolveLocalPaths)
        case CwlOptional(CwlAny) =>
          val (cwlType, cwlValue) = fromIRValue(innerValue,
                                                Some(innerName),
                                                isInput,
                                                fileResolver = fileResolver,
                                                resolveLocalPaths = resolveLocalPaths)
          (CwlOptional.ensureOptional(cwlType), cwlValue)
        case _ => (innerType, inner(innerValue, innerType, innerName))
      }
    }
    innerMulti(value, cwlType, name)
  }

  def fromIR(
      values: Map[DxName, (Type, Value)],
      typeAliases: Map[String, CwlSchema] = Map.empty,
      isInput: Boolean,
      fileResolver: FileSourceResolver = FileSourceResolver.get
  ): Map[DxName, (CwlType, CwlValue)] = {
    values.map {
      case (dxName, (t, v)) =>
        val cwlType = fromIRType(t, typeAliases, isInput)
        dxName -> fromIRValueWithType(v, cwlType, dxName.decoded, isInput, fileResolver)
    }
  }

  def prettyFormatType(cwlType: CwlType): String = {
    cwlType match {
      case CwlOptional(t)  => s"${prettyFormatType(t)}?"
      case CwlMulti(types) => s"(${types.map(prettyFormatType).mkString("|")})"
      case CwlAny          => "any"
      case CwlNull         => "null"
      case CwlBoolean      => "boolean"
      case CwlInt          => "int"
      case CwlLong         => "long"
      case CwlFloat        => "float"
      case CwlDouble       => "double"
      case CwlString       => "string"
      case CwlFile         => "File"
      case CwlDirectory    => "Directory"
      case a: CwlArray if a.hasName =>
        a.name
      case a: CwlArray =>
        s"array<${prettyFormatType(a.itemType)}>"
      case r: CwlRecord if r.hasName =>
        r.name
      case r: CwlRecord =>
        s"record<${r.fields.values.map(f => s"${prettyFormatType(f.cwlType)} ${f.name}").mkString(", ")}>"
      case e: CwlEnum if e.hasName =>
        e.name
      case e: CwlEnum =>
        s"enum<${e.symbolNames.mkString(",")}>"
    }
  }

  def prettyFormatPrimitiveValue(value: PrimitiveValue,
                                 verbose: Boolean = false,
                                 withType: Boolean = false): String = {
    val pretty = value match {
      case BooleanValue(b)              => b.toString
      case IntValue(i)                  => i.toString
      case LongValue(l)                 => l.toString
      case FloatValue(f)                => f.toString
      case DoubleValue(d)               => d.toString
      case StringValue(s)               => s
      case f: FileValue if verbose      => f.toJson.prettyPrint
      case f: FileValue                 => f.toString
      case d: DirectoryValue if verbose => d.toJson.prettyPrint
      case d: DirectoryValue            => d.toString
      case _ =>
        throw new Exception(s"unexpected primitive value ${value}")
    }
    if (withType) {
      s"${prettyFormatType(value.cwlType)}(${pretty})"
    } else {
      pretty
    }
  }

  def prettyFormatValue(value: CwlValue,
                        verbose: Boolean = false,
                        indent: Int = 0,
                        withType: Boolean = false): String = {
    val indentStr = " " * indent
    val pretty = value match {
      case NullValue         => "null"
      case p: PrimitiveValue => prettyFormatPrimitiveValue(p, verbose, withType)
      case ArrayValue(items) =>
        val pretty = if (items.isEmpty) {
          "[]"
        } else if (verbose) {
          s"[\n${indentStr}${items
            .map(prettyFormatValue(_, verbose = true, indent + 2))
            .mkString(s"\n")}\n${indentStr}]"
        } else {
          s"[${items.map(prettyFormatValue(_))}]"
        }
        if (withType) {
          s"array(${pretty})"
        } else {
          pretty
        }
      case ObjectValue(fields) =>
        val pretty = if (fields.isEmpty) {
          "{}"
        } else if (verbose) {
          val fieldStrs = fields.map {
            case (name, value) =>
              s"${name}: ${prettyFormatValue(value, verbose = true, indent + 2)}"
          }
          s"{\n${indentStr}${fieldStrs.mkString("\n")}\n${indentStr}}"
        } else {
          val fieldStrs = fields.map {
            case (name, value) => s"${name}: ${prettyFormatValue(value)}"
          }
          s"{${fieldStrs.mkString(",")}}"
        }
        if (withType) {
          s"object(${pretty})"
        } else {
          pretty
        }
      case _ => throw new Exception(s"unrecognized CWL value ${value})")
    }
    s"${indentStr}${pretty}"
  }

  def prettyFormatEnv(env: Map[String, (CwlType, CwlValue)],
                      verbose: Boolean = false,
                      indent: Int = 2): String = {
    env
      .map {
        case (name, (t, v)) =>
          s"${" " * indent}${name}: ${prettyFormatType(t)} ${prettyFormatValue(v, verbose)}"
      }
      .mkString("\n")
  }

  /**
    * Does a WorkflowStep represent a simple call - i.e. with no
    * scatter or conditional?
    * @param step the WorkflowStep
    * @return
    */
  def isSimpleCall(step: WorkflowStep): Boolean = {
    step.scatter.isEmpty && step.when.isEmpty && step.inputs.forall { inp =>
      // * if there is a default value and a source, they can't both go in the value
      // for the stage input, so a fragment is required
      // * if there is more than one source, a fragment is required to merge them
      // * if there is a valueFrom, a fragment is required to evaluate it
      if (inp.default.nonEmpty) {
        inp.sources.isEmpty
      } else {
        inp.sources.size <= 1
      } &&
      inp.linkMerge.isEmpty &&
      inp.pickValue.isEmpty &&
      inp.valueFrom.isEmpty
    }
  }

  def isArray(t: CwlType): Boolean = {
    t match {
      case _: CwlArray => true
      case CwlMulti(types) =>
        types.exists {
          case _: CwlArray => true
          case _           => false
        }
      case _ => false
    }
  }

  /**
    * Returns true if a down-cast is required for `from` to be compatible
    * with `to`. A down-cast is required to go from a wider type to a more
    * narrow type, for example `Any` -> `File`. Throws an exception of the
    * types are not compatible at all.
    */
  def requiresDowncast(from: CwlType, to: CwlType): Boolean = {
    (from, to) match {
      case (CwlAny, CwlAny)                                                           => false
      case (CwlAny, _)                                                                => true
      case (_: CwlMulti, CwlAny)                                                      => false
      case (fromMulti: CwlMulti, toMulti: CwlMulti) if fromMulti.coercibleTo(toMulti) => false
      case (_: CwlMulti, _: CwlMulti) =>
        throw new Exception(s"type ${from} cannot be downcast to ${to}")
      case (fromMulti: CwlMulti, _) if fromMulti.coercibleTo(to) => true
      case (_: CwlMulti, _) =>
        throw new Exception(s"type ${from} cannot be downcast to ${to}")
      case _ if from.coercibleTo(to) => false
      case _ =>
        throw new Exception(s"type ${from} cannot be downcast to ${to}")
    }
  }

  def isDxPath(path: PathValue): Boolean = {
    path.location.orElse(path.path).exists(_.startsWith(DxPath.DxUriPrefix))
  }

  def toJson(values: Map[DxName, (CwlType, CwlValue)]): JsObject = {
    JsObject(values.map {
      case (dxName, (_, v)) => dxName.decoded -> v.toJson
    })
  }

  /**
    * Format an identifier as the target workflow step to execute. A target is in the format
    * {proc}#{path}, where proc is the top-level process name and path is the path from the
    * top-level process to the target. For example, if the top-level workflow is 'wf' and we want
    * to run a step in a sub-workflow, it might be 'wf#step1/nested/step2'.
    * @param id the identifier to format
    * @param parent the process parent - if non-empty, is prepended to the identifier before
    *               splitting off the top-level workflow name.
    * @return
    */
  def formatTarget(id: Identifier, parent: Option[String] = None): String = {
    val dxName = CwlDxName.fromDecodedName(id.frag)
    val fullDxName = if (parent.nonEmpty) {
      val parentName = dxName.pushDecodedNamespace(parent.get)
      dxName.pushDecodedNamespaces(parentName.getDecodedParts)
    } else {
      dxName
    }
    if (fullDxName.numParts == 1) {
      id.name
    } else {
      val (process, step) = fullDxName.popDecodedNamespace()
      s"${process}#${step.toString}"
    }
  }

  def simplifyProcess(process: Process): Process = {
    process.copySimplifyIds(dropNamespace = true,
                            replacePrefix = (Left(true), None),
                            simplifyAutoNames = false,
                            dropCwlExtension = false)
  }

  def createRuntime(workerPaths: DxWorkerPaths): Runtime = {
    Runtime.create(
        outdir = workerPaths.getOutputFilesDir(ensureExists = true).asJavaPath,
        tmpdir = workerPaths.getTempDir(ensureExists = true).asJavaPath
    )
  }

  def createEvaluatorContext(selfValue: CwlValue = NullValue,
                             selfType: Option[CwlType] = None,
                             selfParam: Option[Identifiable with Loadable] = None,
                             env: Map[String, (CwlType, CwlValue)] = Map.empty,
                             inputParameters: Map[String, InputParameter] = Map.empty,
                             inputDir: Path = Paths.get("."),
                             fileResolver: FileSourceResolver = FileSourceResolver.get,
                             runtime: Runtime = Runtime.empty): EvaluatorContext = {
    val selfFinal = selfParam
      .map { param =>
        val t = selfType.getOrElse {
          param match {
            case p: Parameter => p.cwlType
            case _ =>
              throw new Exception(s"cannot determine self type from paramemter ${param}")
          }
        }
        EvaluatorContext.finalizeInputValue(selfValue, t, param, inputDir, fileResolver)
      }
      .getOrElse(selfValue)
    val valuesFinal = env
      .map {
        case (key, (_, value)) if inputParameters.contains(key) =>
          val param = inputParameters(key)
          key -> EvaluatorContext
            .finalizeInputValue(value, param.cwlType, param, inputDir, fileResolver)
        case (key, (_, value)) => key -> value
      }
      .to(SeqMap)
    EvaluatorContext(selfFinal, ObjectValue(valuesFinal), runtime)
  }
}
