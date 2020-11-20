package dx.core.languages.wdl

import wdlTools.eval.{EvalException, FunctionContext, UserDefinedFunctionImplFactory, WdlValues}
import wdlTools.types.{TypeUtils, UserDefinedFunctionPrototype, WdlTypes}

import scala.collection.immutable.SeqMap

object ArchiveFunction extends UserDefinedFunctionPrototype with UserDefinedFunctionImplFactory {
  val ArchiveFunctionPrefix = "archive_"
  private val deleteSourceFilesType = WdlTypes.T_Optional(WdlTypes.T_Boolean)
  private val outputType = WdlTypes.T_File

  override def getPrototype(funcName: String,
                            inputTypes: Vector[WdlTypes.T]): Option[WdlTypes.T_Function] = {
    if (!funcName.startsWith(ArchiveFunctionPrefix)) {
      None
    } else {
      inputTypes match {
        case Vector(t) if !TypeUtils.isOptional(t) =>
          Some(WdlTypes.T_Function1(funcName, t, outputType))
        case Vector(t, deleteSourceFilesType) =>
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
          Some((WdlTypes.T_Function1(taskName, t, outputType), inputNames))
        case Vector((t, false), (WdlTypes.T_Boolean, true)) =>
          Some((WdlTypes.T_Function2(taskName, t, deleteSourceFilesType, outputType), inputNames))
        case _ => None
      }
    } else {
      None
    }
  }

  def createArchive(ctx: FunctionContext): WdlValues.V = {
    val removeSourceFiles = ctx.args match {
      case Vector(_)                         => false
      case Vector(_, WdlValues.V_Boolean(b)) => b
      case _ =>
        throw new EvalException(s"invalid arguments to archive_* function ${ctx.args}", ctx.loc)
    }

  }

  override def getImpl(funcName: String,
                       args: Vector[WdlValues.V]): Option[FunctionContext => WdlValues.V] = {
    if (funcName.startsWith("archive_")) {
      args match {
        case Vector(_)                         => Some(createArchive)
        case Vector(_, WdlValues.V_Boolean(_)) => Some(createArchive)
        case _                                 => None
      }
    } else {
      None
    }
  }
}
