package dx.core.languages.cwl

import dx.core.ir.Value
import dx.core.ir.Value._
import dx.cwl._

import scala.annotation.tailrec

object CwlUtils {
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

  def fromIR(values: Map[String, Value]): Map[String, (CwlType, CwlValue)] = {
    values.map {
      case (name, value) => name -> fromIRValue(value, Some(name))
    }
  }
}
