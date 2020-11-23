package dx.core.languages.wdl

import java.nio.file.Paths

import dx.core.ir.{LocalizedArchive, PackedArchive}
import wdlTools.eval.{EvalException, FunctionContext, UserDefinedFunctionImplFactory, WdlValues}
import wdlTools.syntax.SourceLocation
import wdlTools.types.{TypeUtils, UserDefinedFunctionPrototype, WdlTypes}

import scala.collection.immutable.SeqMap
import scala.util.matching.Regex

abstract class GenericUserDefinedFunction(pattern: Regex)
    extends UserDefinedFunctionPrototype
    with UserDefinedFunctionImplFactory {
  private var prototypeCache: Map[String, Option[WdlTypes.T_Function]] = Map.empty
  private var taskProxyCache: Map[String, Option[(WdlTypes.T_Function, Vector[String])]] = Map.empty

  def nameMatches(name: String): Boolean = pattern.matches(name)

  protected def createPrototype(funcName: String,
                                inputTypes: Vector[WdlTypes.T]): Option[WdlTypes.T_Function]

  override def getPrototype(funcName: String,
                            inputTypes: Vector[WdlTypes.T]): Option[WdlTypes.T_Function] = {
    if (nameMatches(funcName)) {
      if (prototypeCache.contains(funcName)) {
        prototypeCache(funcName)
      } else {
        val prototype = createPrototype(funcName, inputTypes)
        prototypeCache += (funcName -> prototype)
        prototype
      }
    } else {
      None
    }
  }

  protected def createTaskProxyFunction(
      taskName: String,
      input: SeqMap[String, (WdlTypes.T, Boolean)],
      output: SeqMap[String, WdlTypes.T]
  ): Option[(WdlTypes.T_Function, Vector[String])]

  override def getTaskProxyFunction(
      taskName: String,
      input: SeqMap[String, (WdlTypes.T, Boolean)],
      output: SeqMap[String, WdlTypes.T]
  ): Option[(WdlTypes.T_Function, Vector[String])] = {
    if (nameMatches(taskName)) {
      if (taskProxyCache.contains(taskName)) {
        taskProxyCache(taskName)
      } else {
        val taskProxy = createTaskProxyFunction(taskName, input, output)
        taskProxyCache += (taskName -> taskProxy)
        taskProxy
      }
    } else {
      None
    }
  }

  protected def createImpl(prototype: WdlTypes.T_Function,
                           args: Vector[WdlValues.V],
                           loc: SourceLocation): Option[FunctionContext => WdlValues.V]

  override def getImpl(funcName: String,
                       args: Vector[WdlValues.V],
                       loc: SourceLocation): Option[FunctionContext => WdlValues.V] = {
    if (nameMatches(funcName)) {
      val prototype = prototypeCache
        .get(funcName)
        .flatten
        .orElse(taskProxyCache.get(funcName).flatten.map(_._1))
        .getOrElse(
            throw new EvalException(s"no prototype for generic function ${funcName}")
        )
      createImpl(prototype, args, loc)
    } else {
      None
    }
  }
}

object ArchiveFunction extends GenericUserDefinedFunction("archive_.*".r) {
  private val deleteSourceFilesType = WdlTypes.T_Optional(WdlTypes.T_Boolean)
  private val outputType = WdlTypes.T_File

  override protected def createPrototype(
      funcName: String,
      inputTypes: Vector[WdlTypes.T]
  ): Option[WdlTypes.T_Function] = {
    inputTypes match {
      case Vector(t) if !TypeUtils.isOptional(t) =>
        Some(WdlTypes.T_Function1(funcName, t, outputType))
      case Vector(t, deleteSourceFilesType) =>
        Some(WdlTypes.T_Function2(funcName, t, deleteSourceFilesType, outputType))
      case _ => None
    }
  }

  override protected def createTaskProxyFunction(
      taskName: String,
      input: SeqMap[String, (WdlTypes.T, Boolean)],
      output: SeqMap[String, WdlTypes.T]
  ): Option[(WdlTypes.T_Function, Vector[String])] = {
    if (output.values.toVector == Vector(WdlTypes.T_File)) {
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

  override protected def createImpl(
      prototype: WdlTypes.T_Function,
      args: Vector[WdlValues.V],
      loc: SourceLocation
  ): Option[FunctionContext => WdlValues.V] = {
    val inputType = prototype match {
      case f: WdlTypes.T_Function1 => f.input
      case f: WdlTypes.T_Function2 => f.arg1
      case _ =>
        throw new EvalException(s"invalid prototype for generic archive function ${prototype.name}",
                                loc)
    }
    args match {
      case Vector(_) =>
        Some(PackArchive(prototype.name, inputType).apply)
      case Vector(_, WdlValues.V_Boolean(_)) =>
        Some(PackArchive(prototype.name, inputType).apply)
      case _ => None
    }
  }
}

object UnarchiveFunction extends GenericUserDefinedFunction("unarchive_.*".r) {
  private val inputType = WdlTypes.T_File

  override protected def createPrototype(
      funcName: String,
      inputTypes: Vector[WdlTypes.T]
  ): Option[WdlTypes.T_Function] = {
    if (inputTypes.size == 1 && inputTypes.head == inputType) {
      Some(WdlTypes.T_Function1(funcName, inputType, WdlTypes.T_Any))
    } else {
      None
    }
  }

  override protected def createTaskProxyFunction(
      taskName: String,
      input: SeqMap[String, (WdlTypes.T, Boolean)],
      output: SeqMap[String, WdlTypes.T]
  ): Option[(WdlTypes.T_Function, Vector[String])] = {
    if (input.size == 1 &&
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

  override protected def createImpl(
      prototype: WdlTypes.T_Function,
      args: Vector[WdlValues.V],
      loc: SourceLocation
  ): Option[FunctionContext => WdlValues.V] = {
    val inputType = prototype match {
      case p: WdlTypes.T_Function1 => p.input
      case _ =>
        throw new EvalException(
            s"invalid prototype for generic unarchive function ${prototype.name}",
            loc
        )
    }
    (inputType, args) match {
      case (WdlTypes.T_File, Vector(WdlValues.V_File(_))) => Some(unpackArchive)
      case _                                              => None
    }
  }
}
