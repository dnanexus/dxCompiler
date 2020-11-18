package dx.core.languages.cwl

import java.nio.file.Path

import dx.cwl.{ArrayValue, BooleanValue, CwlArray, CwlBoolean, CwlDirectory, CwlDouble, CwlEnum, CwlFile, CwlFloat, CwlInt, CwlLong, CwlNull, CwlOptional, CwlRecord, CwlString, CwlType, CwlValue, DirectoryValue, DoubleValue, FileValue, FloatValue, IntValue, LongValue, NullValue, Parser, Process, StringValue}
import dx.core.ir.{Type, Value}
import dx.core.ir.Type.{TArray, TBoolean, TDirectory, TFile, TFloat, THash, TInt, TOptional, TSchema, TString}
import dx.core.ir.Value.{VArray, VBoolean, VDirectory, VFile, VFloat, VInt, VNull, VString}
import wdlTools.syntax.SourceLocation
import wdlTools.util.Logger

object CwlUtils {
  val locPlaceholder: SourceLocation = SourceLocation.empty

  def parseSource(sourceCode: Path,
                  logger: Logger = Logger.get): Process = {
    try {
      Parser.parse(sourceCode)
    } catch {
      case e: Exception => {
        logger.error(s"Code is syntactically invalid. (${sourceCode})")
        throw e
      }
    }
  }

  //  def serializeType(t: CwlType): JsValue = {
  //    t match {
  //      case CwlBoolean => JsString("Boolean")
  //      case CwlInt => JsString("Int")
  //      case CwlLong => JsString("Int")
  //      case CwlFloat => JsString("Float")
  //      case CwlDouble => JsString("Float")
  //      case CwlString => JsString("String")
  //      case CwlFile => JsString("File")
  //      case CwlDirectory => JsString("Directory")
  //      case CwlNull => JsString("Null")
  //      // TODO: add T_Array, T_Enum and T_record
  //      case _ =>
  //        throw new Exception(s"Unhandled type ${t}")
  //    }
  //  }
  //  def simpleFromString(s: String): CwlType = {
  //    // recordfield,enum, record, array
  //    s match {
  //      case "Boolean" => CwlBoolean
  //      case "Int" => CwlInt
  //      case "Long" => CwlLong
  //      case "Float" => CwlFloat
  //      case "Double" => CwlDouble
  //      case "String" => CwlString
  //      case "File" => CwlFile
  //      case "Directory" => CwlDirectory
  //      case "null" => CwlNull
  //      case s if s.contains("[") =>
  //        throw new Exception(s"type ${s} is not primitive")
  //      case _ =>
  //        throw UnknownTypeException(s"Unknown type ${s}")
  //    }
  //  }
  //

  def createTSchemaFromRecord(cwlRecord: CwlRecord): TSchema = {
    val name = cwlRecord.name.getOrElse("unknown") // FIXME can name really be empty? if so, which name do we use?
    val fields = cwlRecord.fields.map({ case (key, value) => key -> toIRType(value.types) }) // FIXME: RecordField has name as well
    TSchema(name, fields)
  }

  def createTSchemaFromEnum(cwlEnum: CwlEnum): TSchema = {
    val name = cwlEnum.name.get
    val x = List.fill(cwlEnum.symbols.size)(TString)
    TSchema(name, cwlEnum.symbols.zip(x).toMap)
  }


  // TODO: start only with simple types, then add complex types and then add multiple types
  def toIRType(cwlTypes: Vector[CwlType]): Type = {
    cwlTypes match {
      case _ if cwlTypes.size > 1 && cwlTypes.contains(CwlNull) => TOptional(toIRType(cwlTypes.filter(_ != CwlNull)))
      case _ if cwlTypes.length == 1 => cwlTypes.head match {
        //      case CwlAny => ? // FIXME: this type is missing in WDL even though T_Any exists
        case CwlBoolean => TBoolean
        case CwlInt => TInt
        case CwlLong => TInt
        case CwlDouble => TFloat
        case CwlFloat => TFloat
        case CwlString => TString
        case CwlFile => TFile
        case CwlDirectory => TDirectory
        case CwlNull => TSchema(name = "null", members = Map.empty[String, Type])
        case CwlArray(itemTypes, _, _, _, _) => TArray(toIRType(itemTypes))
        case r: CwlRecord => createTSchemaFromRecord(r)
        case e: CwlEnum => createTSchemaFromEnum(e)
        case _ => throw new Exception(s"Cannot convert CWL type ${cwlTypes.head} to IR")
      }
      case other => throw new NotImplementedError(s"Multiple types are not supported yet! (${other} was given)")
    }
  }


  def toIRTypeMap(cwlTypes: Map[String, Vector[CwlType]]): Map[String, Type] = {
    cwlTypes.map {
      case (name, cwlType) => name -> toIRType(cwlType)
    }
  }

  def fromIRType(irType: Type): CwlType = {
    def inner(innerType: Type): CwlType = {
      innerType match {
        case TBoolean => CwlBoolean
        case TInt => CwlInt
        case TFloat => CwlFloat
        case TString => CwlString
        case TFile => CwlFile
        case TDirectory => CwlDirectory
        case TOptional(t) => inner(t) // TODO:: Needs to return CwlTypes - Vector
        case TSchema(name, _) =>
          throw new Exception(s"Unknown type ${name}")
        case _ =>
          throw new Exception(s"Cannot convert IR type ${innerType} to WDL")
      }
    }

    inner(irType)
  }


  def getDefaultIRValue(cwlTypes: Vector[CwlType]): Value = {
    getDefaultIRValue(toIRType(cwlTypes))
  }

  def getDefaultIRValue(IRType: Type): Value = {
    IRType match {
      case TInt => VInt(0)
      case TFloat => VFloat(0.0)
      case TBoolean => VBoolean(true)
      case TString => VString("")
      case TFile => VFile("placeholder.txt")
      case TOptional(_) => VNull // FIXME ??
      case TArray(_, _) => VArray(Vector.empty[Value])
      case TDirectory => VDirectory(".") // FIXME ??
      case THash => throw new NotImplementedError("THash is not implemented yet")
      case TSchema(_, _) => throw new NotImplementedError("TSchema is not implemented yet")
      case other => throw new NotImplementedError(s"${other} is not supported.")

      // FIXME: add THASH and TSCHEMA
    }
  }

  def getDefaultCWLValue(cwlTypes: Vector[CwlType]): CwlValue = {
    val cwlType = cwlTypes.size match {
      case 1 => cwlTypes.head
      case 0 => throw new Exception("Variable in CWL has to have at least one possible type.")
      case _ => cwlTypes.find(_ != CwlNull).get // Cannot be empty, as there has to be at least 2 types and both cannot be CwlNull
    }
    cwlType match {
      case CwlInt => IntValue(0)
      case CwlFloat => FloatValue(0.0)
      case CwlDouble => DoubleValue(0.0)
      case CwlBoolean => BooleanValue(true)
      case CwlString => StringValue("")
      case CwlFile => FileValue("placeholder.txt")
      case CwlDirectory => DirectoryValue(".") // FIXME ??
      case _: CwlOptional => NullValue // FIXME ??
      case _: CwlArray => ArrayValue(Vector.empty[CwlValue])
      case c: CwlEnum => StringValue(c.symbols.head) // FIXME ??
      case _: CwlRecord => NullValue // FIXME ?
      case other => throw new NotImplementedError(s"${other} is not supported.")
    }
  }

  def toIRValue(cwlValue: CwlValue): Value = {
    cwlValue match {
      case NullValue => VNull
      case BooleanValue(b) => VBoolean(value = b)
      case IntValue(i) => VInt(i)
      case LongValue(l) => VInt(l)
      case FloatValue(f) => VFloat(f)
      case DoubleValue(d) => VFloat(d)
      case StringValue(s) => VString(s)
      case ArrayValue(a) => VArray(a.map(toIRValue))
      case fileValue: FileValue => VFile(fileValue.location.get) // FIXME: default value for location?
      case null => VNull
      case _ => throw new Exception(s"Invalid CWL value ${cwlValue})")
    }
  }


  def toIR(cwl: Map[String, (Vector[CwlType], Option[CwlValue])]): Map[String, (Type, Value)] = {
    cwl.map {
      case (name, (cwlTypes, cwlValue)) =>
        val irType = toIRType(cwlTypes)
        val irValue = cwlValue match {
          case Some(v) => toIRValue(v)
          case None => null
        }
        name -> (irType, irValue)
    }
  }
}


