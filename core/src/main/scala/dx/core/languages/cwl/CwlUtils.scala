package dx.core.languages.cwl

import dx.core.io.DxWorkerPaths
import dx.core.ir.{Type, Value}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.cwl._

import scala.annotation.tailrec

object CwlUtils {
  def toIRSchema(cwlRecord: CwlRecord): TSchema = {
    if (cwlRecord.name.isEmpty) {
      throw new Exception(s"cannot convert schema without name ${cwlRecord}")
    }
    TSchema(
        cwlRecord.name.get,
        cwlRecord.fields.map {
          case (key, value) if value.types.size == 1 =>
            key -> toIRType(value.types.head)
          case _ =>
            throw new Exception("Multi-type fields are not supported")
        }
    )
  }

  def toIRType(cwlType: CwlType): Type = {
    cwlType match {
      case CwlBoolean   => TBoolean
      case CwlInt       => TInt
      case CwlLong      => TInt
      case CwlDouble    => TFloat
      case CwlFloat     => TFloat
      case CwlString    => TString
      case CwlFile      => TFile
      case CwlDirectory => TDirectory
      case a: CwlArray if a.itemTypes.size == 1 =>
        TArray(toIRType(a.itemTypes.head))
      case r: CwlRecord if r.name.isDefined =>
        toIRSchema(r)
      case e: CwlEnum => TEnum(e.symbols)
      case _ =>
        throw new Exception(s"Cannot convert CWL type ${cwlType} to IR")
    }
  }

  def toIRValue(cwlValue: CwlValue): Value = {
    cwlValue match {
      case NullValue         => VNull
      case BooleanValue(b)   => VBoolean(value = b)
      case IntValue(i)       => VInt(i)
      case LongValue(l)      => VInt(l)
      case FloatValue(f)     => VFloat(f)
      case DoubleValue(d)    => VFloat(d)
      case StringValue(s)    => VString(s)
      case f: FileValue      => VFile(f.toString)
      case d: DirectoryValue => VDirectory(d.toString)
      case ArrayValue(a)     => VArray(a.map(toIRValue))
      case ObjectValue(m) =>
        VHash(m.map {
          case (key, value) => key -> toIRValue(value)
        })
      case _ => throw new Exception(s"Invalid CWL value ${cwlValue})")
    }
  }

  def toIRValue(cwlValue: CwlValue, cwlType: CwlType): Value = {
    (cwlType, cwlValue) match {
      case (CwlOptional(_), NullValue)       => VNull
      case (CwlOptional(t), _)               => toIRValue(cwlValue, t)
      case (CwlBoolean, BooleanValue(b))     => VBoolean(b)
      case (CwlInt, IntValue(i))             => VInt(i)
      case (CwlLong, LongValue(l))           => VInt(l)
      case (CwlFloat, FloatValue(f))         => VFloat(f)
      case (CwlDouble, DoubleValue(d))       => VFloat(d)
      case (t: CwlNumber, n: NumericValue)   => toIRValue(n.coerceTo(t), t)
      case (CwlString, StringValue(s))       => VString(s)
      case (CwlFile, f: FileValue)           => VFile(f.toString)
      case (CwlFile, StringValue(s))         => VFile(s)
      case (CwlDirectory, d: DirectoryValue) => VDirectory(d.toString)
      case (CwlDirectory, StringValue(s))    => VDirectory(s)
      case (array: CwlArray, ArrayValue(items)) if array.itemTypes.size == 1 =>
        VArray(items.map(toIRValue(_, array.itemTypes.head)))
      case (record: CwlRecord, ObjectValue(members)) =>
        VHash(members.map {
          case (name, value)
              if record.fields.contains(name) && record.fields(name).types.size == 1 =>
            name -> toIRValue(value, record.fields(name).types.head)
        })
      case (enum: CwlEnum, StringValue(s)) if enum.symbols.contains(s) =>
        VString(s)
      case _ => throw new Exception(s"Invalid CWL value ${cwlValue})")
    }
  }

  def toIR(cwl: Map[String, (CwlType, CwlValue)]): Map[String, (Type, Value)] = {
    cwl.map {
      case (name, (cwlType, cwlValue)) =>
        val irType = toIRType(cwlType)
        val irValue = toIRValue(cwlValue, cwlType)
        name -> (irType, irValue)
    }
  }

  def fromIRType(irType: Type, typeAliases: Map[String, CwlSchema] = Map.empty): CwlType = {
    irType match {
      case TOptional(t) => CwlOptional(fromIRType(t))
      case TBoolean     => CwlBoolean
      case TInt         => CwlLong
      case TFloat       => CwlDouble
      case TString      => CwlString
      case TFile        => CwlFile
      case TDirectory   => CwlDirectory
      case TArray(t, _) => CwlArray(Vector(fromIRType(t, typeAliases)))
      case TSchema(name, _) if typeAliases.contains(name) =>
        typeAliases(name)
      case TSchema(name, members) =>
        CwlRecord(members.map {
          case (name, t) => name -> CwlRecordField(name, Vector(fromIRType(t, typeAliases)))
        }, Some(name))
      case TEnum(allowedValues) => CwlEnum(allowedValues)
      case _ =>
        throw new Exception(s"Cannot convert IR type ${irType} to CWL")
    }
  }

  def fromIRValue(value: Value, name: Option[String]): (CwlType, CwlValue) = {
    value match {
      case VNull         => (CwlNull, NullValue)
      case VBoolean(b)   => (CwlBoolean, BooleanValue(b))
      case VInt(i)       => (CwlLong, LongValue(i))
      case VFloat(f)     => (CwlDouble, DoubleValue(f))
      case VString(s)    => (CwlString, StringValue(s))
      case VFile(f)      => (CwlFile, FileValue(f))
      case VDirectory(d) => (CwlDirectory, DirectoryValue(d))
      case VArray(array) =>
        val (types, values) = array.zipWithIndex.map {
          case (v, i) => fromIRValue(v, name.map(n => s"${n}[${i}]"))
        }.unzip
        (CwlArray(types.distinct), ArrayValue(values))
      case VHash(fields) =>
        val (types, values) = fields.map {
          case (key, value) =>
            val (cwlType, cwlValue) = fromIRValue(value, name.map(n => s"${n}[${key}]"))
            (key -> cwlType, key -> cwlValue)
        }.unzip
        // create an anonymous record schema
        val schemaType = CwlRecord(types.map {
          case (name, t) => name -> CwlRecordField(name, types = Vector(t))
        }.toMap)
        (schemaType, ObjectValue(values.toMap))
      case _ =>
        throw new Exception(
            s"cannot convert ${name.getOrElse("IR")} value ${value} to WDL value"
        )
    }
  }

  def fromIR(values: Map[String, Value]): Map[String, (CwlType, CwlValue)] = {
    values.map {
      case (name, value) => name -> fromIRValue(value, Some(name))
    }
  }

  def fromIRValue(value: Value, cwlTypes: Vector[CwlType], name: String): (CwlType, CwlValue) = {
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
            case (item, i) => fromIRValue(item, array.itemTypes, s"${innerName}[${i}]")._2
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
          ObjectValue(members.map {
            case (key, value) =>
              key -> fromIRValue(value, record.fields(key).types, s"${innerName}[${key}]")._2
          })
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
          fromIRValue(value, Some(name))
        } else {
          throw new Exception(s"Cannot convert ${name} (${cwlTypes}, ${value}) to CWL value")
        }
      }
  }

  def fromIR(values: Map[String, (Type, Value)]): Map[String, (CwlType, CwlValue)] = {
    values.map {
      case (name, (t, v)) =>
        val cwlType = fromIRType(t)
        name -> fromIRValue(v, Vector(cwlType), name)
    }
  }

  def createRuntime(workerPaths: DxWorkerPaths): Runtime = {
    Runtime.create(
        outdir = workerPaths.getOutputFilesDir(ensureExists = true),
        tmpdir = workerPaths.getTempDir(ensureExists = true)
    )
  }

  def createEvauatorContext(runtime: Runtime,
                            env: Map[String, (CwlType, CwlValue)] = Map.empty,
                            self: CwlValue = NullValue): EvaluatorContext = {
    val values = env.map {
      case (key, (_, value)) => key -> value
    }
    EvaluatorContext(self, ObjectValue(values), runtime)
  }
}
