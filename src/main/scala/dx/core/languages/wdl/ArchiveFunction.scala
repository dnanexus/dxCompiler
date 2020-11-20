package dx.core.languages.wdl

import java.nio.file.Paths

import dx.core.ir.{LocalizedArchive, PackedArchive}
import wdlTools.eval.{EvalException, FunctionContext, UserDefinedFunctionImplFactory, WdlValues}
import wdlTools.types.{TypeUtils, UserDefinedFunctionPrototype, WdlTypes}

import scala.collection.immutable.SeqMap

case class ArchiveFunction(initialInputValueTypes: Map[String, WdlTypes.T] = Map.empty)
    extends UserDefinedFunctionPrototype
    with UserDefinedFunctionImplFactory {
  val ArchiveFunctionPrefix = "archive_"
  private val deleteSourceFilesType = WdlTypes.T_Optional(WdlTypes.T_Boolean)
  private val outputType = WdlTypes.T_File
  private var inputValueTypes = initialInputValueTypes

  override def getPrototype(funcName: String,
                            inputTypes: Vector[WdlTypes.T]): Option[WdlTypes.T_Function] = {
    if (!funcName.startsWith(ArchiveFunctionPrefix)) {
      None
    } else {
      inputTypes match {
        case Vector(t) if !TypeUtils.isOptional(t) =>
          inputValueTypes += (funcName -> t)
          Some(WdlTypes.T_Function1(funcName, t, outputType))
        case Vector(t, deleteSourceFilesType) =>
          inputValueTypes += (funcName -> t)
          Some(WdlTypes.T_Function2(funcName, t, deleteSourceFilesType, outputType))
        case _ => None
      }
    }
  }

  override def getTaskProxyFunction(
      taskName: String,
      input: SeqMap[String, (WdlTypes.T, Boolean)],
      output: SeqMap[String, WdlTypes.T]
  ): Option[(WdlTypes.T_Function, Vector[String])] = {
    if (taskName
          .startsWith(ArchiveFunctionPrefix) && output.values.toVector == Vector(WdlTypes.T_File)) {
      val (inputNames, inputTypes) = input.toVector.unzip
      inputTypes match {
        case Vector((t, false)) =>
          inputValueTypes += (taskName -> t)
          Some((WdlTypes.T_Function1(taskName, t, outputType), inputNames))
        case Vector((t, false), (WdlTypes.T_Boolean, true)) =>
          inputValueTypes += (taskName -> t)
          Some((WdlTypes.T_Function2(taskName, t, deleteSourceFilesType, outputType), inputNames))
        case _ => None
      }
    } else {
      None
    }
  }

  case class PackArchive(funcName: String, inputValueType: WdlTypes.T) {
    def apply(ctx: FunctionContext): WdlValues.V = {
      val removeSourceFiles = ctx.args match {
        case Vector(_)                         => false
        case Vector(_, WdlValues.V_Boolean(b)) => b
        case _ =>
          throw new EvalException(s"invalid arguments to archive_* function ${ctx.args}", ctx.loc)
      }
      val irType = WdlUtils.toIRType(inputValueType)
      val irValue = WdlUtils.toIRValue(ctx.args(0), inputValueType)
      val parentDir = ctx.paths.getRootDir(true)
      val archive = LocalizedArchive(irType, irValue)(parentDir = Some(parentDir))
      val packedArchive = archive.pack(removeSourceFiles)
      WdlValues.V_File(packedArchive.path.toString)
    }
  }

  override def getImpl(funcName: String,
                       args: Vector[WdlValues.V]): Option[FunctionContext => WdlValues.V] = {
    if (funcName.startsWith(ArchiveFunctionPrefix)) {
      val inputValueType = inputValueTypes.getOrElse(
          funcName,
          throw new EvalException(s"no prototype for generic archive function ${funcName}")
      )
      args match {
        case Vector(_) => Some(PackArchive(funcName, inputValueType).apply)
        case Vector(_, WdlValues.V_Boolean(_)) =>
          Some(PackArchive(funcName, inputValueType).apply)
        case _ => None
      }
    } else {
      None
    }
  }
}

object UnarchiveFunction extends UserDefinedFunctionPrototype with UserDefinedFunctionImplFactory {
  val UnarchiveFunctionPrefix = "unarchive_"
  private val inputType = WdlTypes.T_File

  override def getPrototype(funcName: String,
                            inputTypes: Vector[WdlTypes.T]): Option[WdlTypes.T_Function] = {
    if (funcName.startsWith(UnarchiveFunctionPrefix) &&
        inputTypes.size == 1 &&
        inputTypes.head == inputType) {
      Some(WdlTypes.T_Function1(funcName, inputType, WdlTypes.T_Any))
    } else {
      None
    }
  }

  override def getTaskProxyFunction(
      taskName: String,
      input: SeqMap[String, (WdlTypes.T, Boolean)],
      output: SeqMap[String, WdlTypes.T]
  ): Option[(WdlTypes.T_Function, Vector[String])] = {
    if (taskName.startsWith(UnarchiveFunctionPrefix) &&
        input.size == 1 &&
        input.values.head == (inputType, false) &&
        output.size == 1) {
      Some((WdlTypes.T_Function1(taskName, inputType, output.values.head), Vector(input.keys.head)))
    } else {
      None
    }
  }

  def unpackArchive(ctx: FunctionContext): WdlValues.V = {
    val path = ctx.getOneArg match {
      case WdlValues.V_File(path) => Paths.get(path)
      case other =>
        throw new EvalException(s"expected a path argument, not ${other}", ctx.loc)
    }
    val packedArchive = PackedArchive(path)()
    val (localizedArchive, _) = packedArchive.localize()
    val wdlType = WdlUtils.fromIRType(localizedArchive.irType)
    WdlUtils.fromIRValue(localizedArchive.irValue, wdlType, "unknown")
  }

  override def getImpl(funcName: String,
                       args: Vector[WdlValues.V]): Option[FunctionContext => WdlValues.V] = {
    if (funcName.startsWith(UnarchiveFunctionPrefix)) {
      args match {
        case Vector(WdlValues.V_File(_)) => Some(unpackArchive)
        case _                           => None
      }
    } else {
      None
    }
  }
}
