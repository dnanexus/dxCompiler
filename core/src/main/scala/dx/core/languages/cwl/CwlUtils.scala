package dx.core.languages.cwl

import dx.api.DxPath
import dx.core.io.DxWorkerPaths
import dx.core.ir.{Type, Value}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.cwl._
import spray.json._

import scala.annotation.tailrec
import scala.collection.immutable.TreeSeqMap

object CwlUtils {

  /**
    * Attempts to convert a Vector of alternative types to a single type.
    * If the Vector is of size 2 and one of the types is CwlNull, returns
    * the other type, converted to an optional. If the Vector contains
    * CwlAny, then CwlAny is returned. Otherwise throws an Exception
    * unless the Vector contains exactly one type.
    * @param types Vector of CwlTypes
    * @return
    */
  def getSingleType(types: Vector[CwlType]): CwlType = {
    types match {
      case Vector(t)               => t
      case v if v.contains(CwlAny) => CwlAny
      case Vector(CwlNull, other)  => CwlOptional.ensureOptional(other)
      case Vector(other, CwlNull)  => CwlOptional.ensureOptional(other)
      case Vector() =>
        throw new Exception("No type specified")
      case _ =>
        throw new Exception("Multi-type fields are not supported")
    }
  }

  def toIRSchema(cwlRecord: CwlRecord): TSchema = {
    if (cwlRecord.name.isEmpty) {
      throw new Exception(s"cannot convert schema without name ${cwlRecord}")
    }
    TSchema(
        cwlRecord.name.get,
        cwlRecord.fields.map {
          case (key, value) => key -> toIRType(getSingleType(value.types))
        }
    )
  }

  def toIRType(cwlType: CwlType): Type = {
    cwlType match {
      case CwlOptional(t) => TOptional(toIRType(t))
      case CwlBoolean     => TBoolean
      case CwlInt         => TInt
      case CwlLong        => TInt
      case CwlDouble      => TFloat
      case CwlFloat       => TFloat
      case CwlString      => TString
      case CwlFile        => TFile
      case CwlDirectory   => TDirectory
      case a: CwlArray if a.itemTypes.size == 1 =>
        TArray(toIRType(a.itemTypes.head))
      case r: CwlRecord if r.name.isDefined =>
        toIRSchema(r)
      case e: CwlEnum => TEnum(e.symbols)
      case _ =>
        throw new Exception(s"Cannot convert CWL type ${cwlType} to IR")
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
      case f: FileValue      => (TFile, VFile(f.toString))
      case d: DirectoryValue => (TDirectory, VDirectory(d.toString))
      case ArrayValue(items) =>
        val (itemTypes, itemValues, optional) =
          items.foldLeft(Set.empty[Type], Vector.empty[Value], false) {
            case ((types, values, optional), item) =>
              item match {
                case NullValue => (types, values :+ VNull, true)
                case _ =>
                  val (t, v) = toIRValue(item)
                  if (types.isEmpty || types.contains(t)) {
                    (types + t, values :+ v, optional)
                  } else {
                    throw new Exception(s"array ${items} contains values of multiple types")
                  }
              }
          }
        if (itemTypes.isEmpty) {
          throw new Exception("cannot determine type of array with only null items")
        }
        val itemType = itemTypes.head match {
          case t if optional => ensureOptional(t)
          case t             => t
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
               .to(TreeSeqMap)
         ))
      case _ => throw new Exception(s"Invalid CWL value ${cwlValue})")
    }
  }

  def toIRValue(cwlValue: CwlValue, cwlType: CwlType): (Type, Value) = {
    (cwlType, cwlValue) match {
      case (CwlAny, _)                 => toIRValue(cwlValue)
      case (CwlOptional(_), NullValue) => (toIRType(cwlType), VNull)
      case (CwlOptional(t), _) =>
        val (irType, irValue) = toIRValue(cwlValue, t)
        (TOptional(irType), irValue)
      case (CwlBoolean, BooleanValue(b))     => (TBoolean, VBoolean(b))
      case (CwlInt, IntValue(i))             => (TInt, VInt(i))
      case (CwlLong, LongValue(l))           => (TInt, VInt(l))
      case (CwlFloat, FloatValue(f))         => (TFloat, VFloat(f))
      case (CwlDouble, DoubleValue(d))       => (TFloat, VFloat(d))
      case (t: CwlNumber, n: NumericValue)   => toIRValue(n.coerceTo(t), t)
      case (CwlString, StringValue(s))       => (TString, VString(s))
      case (CwlFile, f: FileValue)           => (TFile, VFile(f.toString))
      case (CwlFile, StringValue(s))         => (TFile, VFile(s))
      case (CwlDirectory, d: DirectoryValue) => (TDirectory, VDirectory(d.toString))
      case (CwlDirectory, StringValue(s))    => (TDirectory, VDirectory(s))
      case (array: CwlArray, ArrayValue(items)) =>
        val itemType = getSingleType(array.itemTypes)
        val (irItemTypes, irItems, optional) =
          items.foldLeft(Set.empty[Type], Vector.empty[Value], false) {
            case ((irTypes, irItems, _), NullValue) => (irTypes, irItems :+ VNull, true)
            case ((irTypes, irItems, optional), i) =>
              val (t, v) = toIRValue(i, itemType)
              if (irTypes.nonEmpty && !irTypes.contains(t)) {
                throw new Exception(s"array ${array} contains items of multiple types")
              }
              (irTypes + t, irItems :+ v, optional || Type.isOptional(t))
          }
        val irItemType = irItemTypes.headOption.getOrElse(toIRType(itemType))
        val irType = TArray(if (optional) {
          TOptional(irItemType)
        } else {
          irItemType
        })
        (irType, VArray(irItems))
      case (record: CwlRecord, ObjectValue(fields)) =>
        val (types, values) =
          fields.foldLeft(TreeSeqMap.empty[String, Type], TreeSeqMap.empty[String, Value]) {
            case ((types, values), (name, value)) if record.fields.contains(name) =>
              val cwlType = getSingleType(record.fields(name).types)
              val (irType, irValue) = toIRValue(value, cwlType)
              (types + (name -> irType), values + (name -> irValue))
            case (name, _) =>
              throw new Exception(s"invalid field ${name}")
          }
        val irType = if (record.name.isDefined) {
          TSchema(record.name.get, types)
        } else {
          THash
        }
        (irType, VHash(values))
      case (enum: CwlEnum, StringValue(s)) if enum.symbols.contains(s) =>
        (TEnum(enum.symbols), VString(s))
      case _ => throw new Exception(s"Invalid CWL value ${cwlValue})")
    }
  }

  /**
    * The CWL "Any" type does not have a DNAnexus equivalent, so
    * we need to wrap the value (and optionally also the actual
    * type) in a hash.
    */
  def anyToIRValue(cwlValue: CwlValue, cwlType: Option[CwlType] = None): (Type, VHash) = {
    val (irType, irValue) = if (cwlType.isDefined) {
      toIRValue(cwlValue, cwlType.get)
    } else {
      toIRValue(cwlValue)
    }
    (irType, VHash(TreeSeqMap("value" -> irValue)))
  }

  def toIR(cwl: Map[String, (CwlType, CwlValue)]): Map[String, (Type, Value)] = {
    cwl.map {
      case (name, (cwlType, cwlValue)) => name -> toIRValue(cwlValue, cwlType)
    }
  }

  def fromIRType(irType: Type,
                 typeAliases: Map[String, CwlSchema] = Map.empty,
                 isInput: Boolean): CwlType = {
    def inner(innerType: Type): CwlType = {
      innerType match {
        case TOptional(t) => CwlOptional(inner(t))
        case TBoolean     => CwlBoolean
        case TInt         => CwlLong
        case TFloat       => CwlDouble
        case TString      => CwlString
        case TFile        => CwlFile
        case TDirectory   => CwlDirectory
        case TArray(t, _) => CwlArray(Vector(inner(t)))
        case TSchema(name, _) if typeAliases.contains(name) =>
          typeAliases(name)
        case TSchema(name, members) if isInput =>
          CwlInputRecord(members.map {
            case (name, t) => name -> CwlInputRecordField(name, Vector(inner(t)))
          }, Some(name))
        case TSchema(name, members) =>
          CwlOutputRecord(members.map {
            case (name, t) => name -> CwlOutputRecordField(name, Vector(inner(t)))
          }, Some(name))
        case TEnum(allowedValues) => CwlEnum(allowedValues)
        case _ =>
          throw new Exception(s"Cannot convert IR type ${irType} to CWL")
      }
    }
    inner(irType)
  }

  def fromIRValue(value: Value, name: Option[String], isInput: Boolean): (CwlType, CwlValue) = {
    def inner(innerValue: Value, innerName: Option[String]): (CwlType, CwlValue) = {
      innerValue match {
        case VNull         => (CwlNull, NullValue)
        case VBoolean(b)   => (CwlBoolean, BooleanValue(b))
        case VInt(i)       => (CwlLong, LongValue(i))
        case VFloat(f)     => (CwlDouble, DoubleValue(f))
        case VString(s)    => (CwlString, StringValue(s))
        case VFile(f)      => (CwlFile, FileValue(f))
        case VDirectory(d) => (CwlDirectory, DirectoryValue(d))
        case VArray(array) =>
          val (types, values) = array.zipWithIndex.map {
            case (v, i) => inner(v, innerName.map(n => s"${n}[${i}]"))
          }.unzip
          (CwlArray(types.distinct), ArrayValue(values))
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
                    case (name, t) => name -> CwlInputRecordField(name, types = Vector(t))
                  }
                  .to(TreeSeqMap)
            )
          } else {
            CwlOutputRecord(
                types
                  .map {
                    case (name, t) => name -> CwlOutputRecordField(name, types = Vector(t))
                  }
                  .to(TreeSeqMap)
            )
          }
          (schemaType, ObjectValue(values.to(TreeSeqMap)))
        case _ =>
          throw new Exception(
              s"cannot convert ${name.getOrElse("IR")} value ${value} to WDL value"
          )
      }
    }
    inner(value, name)
  }

  def fromIRValues(values: Map[String, Value],
                   isInput: Boolean): Map[String, (CwlType, CwlValue)] = {
    values.map {
      case (name, value) => name -> fromIRValue(value, Some(name), isInput)
    }
  }

  def fromIRValue(value: Value,
                  cwlTypes: Vector[CwlType],
                  name: String,
                  isInput: Boolean): (CwlType, CwlValue) = {
    @tailrec
    def inner(innerValue: Value, innerType: CwlType, innerName: String): CwlValue = {
      (innerType, innerValue) match {
        case (CwlOptional(_) | CwlNull, VNull) => NullValue
        case (CwlOptional(t), _)               => inner(innerValue, t, innerName)
        case (CwlBoolean, VBoolean(b))         => BooleanValue(b)
        case (CwlInt, VInt(i)) if i.isValidInt => IntValue(i)
        case (CwlLong, VInt(l))                => LongValue(l)
        case (CwlFloat, VFloat(f))             => FloatValue(f.toFloat)
        case (CwlFloat, VInt(i))               => FloatValue(i.toFloat)
        case (CwlDouble, VFloat(f))            => FloatValue(f)
        case (CwlDouble, VInt(i))              => FloatValue(i.toDouble)
        case (CwlString, VString(s))           => StringValue(s)
        case (CwlFile, VString(path))          => FileValue(path)
        case (CwlFile, VFile(path))            => FileValue(path)
        case (CwlDirectory, VString(path))     => DirectoryValue(path)
        case (CwlDirectory, VFile(path))       => DirectoryValue(path)
        case (array: CwlArray, VArray(items)) =>
          ArrayValue(items.zipWithIndex.map {
            case (item, i) => fromIRValue(item, array.itemTypes, s"${innerName}[${i}]", isInput)._2
          })
        case (record: CwlRecord, VHash(members)) =>
          // ensure 1) members keys are a subset of memberTypes keys, 2) members
          // values are convertable to the corresponding types, and 3) any keys
          // in memberTypes that do not appear in members are optional
          val keys1 = members.keySet
          val keys2 = record.fields.keySet
          val extra = keys2.diff(keys1)
          if (extra.nonEmpty) {
            throw new Exception(
                s"struct ${record.name} value has members that do not appear in the struct definition: ${extra}"
            )
          }
          val missingNonOptional =
            keys1.diff(keys2).map(key => key -> record.fields(key)).filterNot {
              case (_, field) if CwlOptional.anyOptional(field.types) => false
              case _                                                  => true
            }
          if (missingNonOptional.nonEmpty) {
            throw new Exception(
                s"struct ${record.name} value is missing non-optional members ${missingNonOptional}"
            )
          }
          ObjectValue(
              members.map {
                case (key, value) =>
                  key -> fromIRValue(value,
                                     record.fields(key).types,
                                     s"${innerName}[${key}]",
                                     isInput)._2
              }
          )
        case (enum: CwlEnum, VString(s)) if enum.symbols.contains(s) =>
          StringValue(s)
        case _ =>
          throw new Exception(s"cannot translate ${innerValue} to CwlValue of type ${innerType}")
      }
    }
    cwlTypes.iterator
      .map { t =>
        try {
          Some((t, inner(value, t, name)))
        } catch {
          case _: Throwable => None
        }
      }
      .collectFirst { case Some(result) => result }
      .getOrElse {
        if (cwlTypes.contains(CwlAny)) {
          fromIRValue(value, Some(name), isInput)
        } else {
          throw new Exception(s"Cannot convert ${name} (${cwlTypes}, ${value}) to CWL value")
        }
      }
  }

  def fromIR(values: Map[String, (Type, Value)],
             typeAliases: Map[String, CwlSchema] = Map.empty,
             isInput: Boolean): Map[String, (CwlType, CwlValue)] = {
    values.map {
      case (name, (t, v)) =>
        val cwlType = fromIRType(t, typeAliases, isInput)
        name -> fromIRValue(v, Vector(cwlType), name, isInput)
    }
  }

  def isDxFile(file: FileValue): Boolean = {
    file.location.exists(_.startsWith(DxPath.DxUriPrefix))
  }

  def toJson(values: Map[String, (CwlType, CwlValue)]): JsObject = {
    JsObject(values.map {
      case (name, (_, v)) => name -> v.toJson
    })
  }

  def createRuntime(workerPaths: DxWorkerPaths): Runtime = {
    Runtime.create(
        outdir = workerPaths.getOutputFilesDir(ensureExists = true),
        tmpdir = workerPaths.getTempDir(ensureExists = true)
    )
  }

  def createEvaluatorContext(runtime: Runtime,
                             env: Map[String, (CwlType, CwlValue)] = Map.empty,
                             self: CwlValue = NullValue): EvaluatorContext = {
    val values = env
      .map {
        case (key, (_, value)) => key -> value
      }
      .to(TreeSeqMap)
    EvaluatorContext(self, ObjectValue(values), runtime)
  }
}
