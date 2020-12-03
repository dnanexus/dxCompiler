package dx.core.languages.cwl

import dx.core.ir.Value
import dx.core.ir.Value._
import dx.cwl._

import scala.annotation.tailrec

object CwlUtils {

  def isOptional(t: CwlType): Boolean = {
    t match {
      case CwlOptional(_) => true
      case _              => false
    }
  }

  def anyOptional(types: Vector[CwlType]): Boolean = {
    types.exists(isOptional)
  }
  def fromIRValue(value: Value, name: Option[String]): CwlValue = {
    value match {
      case VNull         => NullValue
      case VBoolean(b)   => BooleanValue(b)
      case VInt(i)       => LongValue(i)
      case VFloat(f)     => DoubleValue(f)
      case VString(s)    => StringValue(s)
      case VFile(f)      => FileValue(f)
      case VDirectory(d) => DirectoryValue(d)
      case VArray(array) =>
        ArrayValue(array.zipWithIndex.map {
          case (v, i) => fromIRValue(v, name.map(n => s"${n}[${i}]"))
        })
      case VHash(fields) =>
        ObjectValue(fields.map {
          case (key, value) => key -> fromIRValue(value, name.map(n => s"${n}[${key}]"))
        })
      case _ =>
        throw new Exception(
            s"cannot convert ${name.getOrElse("IR")} value ${value} to WDL value"
        )
    }
  }

  def fromIRValue(value: Value, cwlTypes: Vector[CwlType], name: String): CwlValue = {
    @tailrec
    def inner(innerValue: Value, innerType: CwlType, innerName: String): Option[CwlValue] = {
      (innerType, innerValue) match {
        case (CwlOptional(_) | CwlNull, VNull) => Some(NullValue)
        case (CwlOptional(t), _)               => inner(innerValue, t, innerName)
        case (CwlBoolean, VBoolean(b))         => Some(BooleanValue(b))
        case (CwlInt, VInt(i)) if i.isValidInt => Some(IntValue(i))
        case (CwlLong, VInt(l))                => Some(LongValue(l))
        case (CwlFloat, VFloat(f))             => Some(FloatValue(f.toFloat))
        case (CwlFloat, VInt(i))               => Some(FloatValue(i.toFloat))
        case (CwlDouble, VFloat(f))            => Some(FloatValue(f))
        case (CwlDouble, VInt(i))              => Some(FloatValue(i.toDouble))
        case (CwlString, VString(s))           => Some(StringValue(s))
        case (CwlFile, VString(path))          => Some(FileValue(path))
        case (CwlFile, VFile(path))            => Some(FileValue(path))
        case (CwlDirectory, VString(path))     => Some(DirectoryValue(path))
        case (CwlDirectory, VFile(path))       => Some(DirectoryValue(path))
        case (array: CwlArray, VArray(items)) =>
          Some(ArrayValue(items.zipWithIndex.map {
            case (item, i) => fromIRValue(item, array.itemTypes, s"${innerName}[${i}]")
          }))
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
              case (_, field) if anyOptional(field.types) => false
              case _                                      => true
            }
          if (missingNonOptional.nonEmpty) {
            throw new Exception(
                s"struct ${record.name} value is missing non-optional members ${missingNonOptional}"
            )
          }
          Some(ObjectValue(members.map {
            case (key, value) =>
              key -> fromIRValue(value, record.fields(key).types, s"${innerName}[${key}]")
          }))
        case (enum: CwlEnum, VString(s)) if enum.symbols.contains(s) =>
          Some(StringValue(s))
        case _ => None
      }
    }
    cwlTypes.iterator
      .map(t => inner(value, t, name))
      .collectFirst {
        case Some(value) => value
      }
      .getOrElse {
        if (cwlTypes.contains(CwlAny)) {
          fromIRValue(value, Some(name))
        } else {
          throw new Exception(s"Cannot convert ${name} (${cwlTypes}, ${value}) to CWL value")
        }
      }
  }
}
