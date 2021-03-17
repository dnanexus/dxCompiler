package dx.executor.wdl

import dx.AppInternalException
import dx.api.{DxExecution, DxFile, DxObject, DxPath, Field}
import dx.core.Constants
import dx.core.ir.{
  Block,
  BlockKind,
  ExecutableLink,
  Manifest,
  Parameter,
  ParameterLink,
  Type,
  Value,
  ValueSerde
}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.core.languages.wdl.{
  ComputedBlockInput,
  OptionalBlockInput,
  OverridableBlockInputWithDynamicDefault,
  OverridableBlockInputWithStaticDefault,
  RequiredBlockInput,
  Runtime,
  VersionSupport,
  WdlBlock,
  WdlBlockInput,
  WdlUtils
}
import dx.executor.{BlockContext, JobMeta, WorkflowExecutor}
import dx.util.{DefaultBindings, FileNode, JsUtils, LocalFileSource, Logger, TraceLevel}
import dx.util.CollectionUtils.IterableOnceExtensions
import dx.util.protocols.DxFileSource
import spray.json._
import wdlTools.eval.{Eval, EvalUtils, Meta, WdlValueBindings}
import wdlTools.eval.WdlValues._
import wdlTools.exec.{InputOutput, TaskInputOutput}
import wdlTools.types.{TypeUtils, TypedAbstractSyntax => TAT}
import wdlTools.types.WdlTypes._

case class BlockIO(block: WdlBlock, logger: Logger)

object WdlWorkflowExecutor {
  def create(jobMeta: JobMeta, separateOutputs: Boolean): WdlWorkflowExecutor = {
    // parse the workflow source code to get the WDL document
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceString(jobMeta.sourceCode, jobMeta.fileResolver)
    val workflow = doc.workflow.getOrElse(
        throw new RuntimeException("This document should have a workflow")
    )
    val tasks = doc.elements.collect {
      case task: TAT.Task => task.name -> task
    }.toMap
    WdlWorkflowExecutor(doc.source,
                        workflow,
                        versionSupport,
                        tasks,
                        typeAliases.toMap,
                        jobMeta,
                        separateOutputs)
  }

  // this method is exposed for unit testing
  def getComplexScatterName(items: Iterator[Option[String]],
                            maxLength: Int = WorkflowExecutor.JobNameLengthLimit): String = {
    // Create a name by concatenating the initial elements of the array.
    // Limit the total size of the name.
    val (_, strings, hasMore) =
      items.foldLeftWhile((-1, Vector.empty[String], false))(_._1 < maxLength) {
        case ((length, strings, _), Some(s)) =>
          val newLength = length + s.length + 1
          if (newLength > maxLength) {
            (newLength, strings, true)
          } else {
            (newLength, strings :+ s, false)
          }
        case ((length, strings, _), None) =>
          (length, strings, false)
      }
    val itemStr = strings.mkString(",")
    if (hasMore) {
      s"${itemStr},..."
    } else {
      itemStr
    }
  }
}

// TODO: implement graph building from workflow - input and output parameters
//  should be ordered in calls to InputOutput.*
/**
  * Executor for WDL workflows.
  * @param docSource the FileNode from which the workflow source originated.
  * @param workflow the Workflow object.
  * @param versionSupport WDL version-specific functions.
  * @param tasks the WDL tasks called by the workflow.
  * @param wdlTypeAliases Struct definitions.
  * @param jobMeta functions to access job and executable metadata.
  */
case class WdlWorkflowExecutor(docSource: FileNode,
                               workflow: TAT.Workflow,
                               versionSupport: VersionSupport,
                               tasks: Map[String, TAT.Task],
                               wdlTypeAliases: Map[String, T_Struct],
                               jobMeta: JobMeta,
                               separateOutputs: Boolean)
    extends WorkflowExecutor[WdlBlock](jobMeta, separateOutputs) {
  private val logger = jobMeta.logger
  private lazy val evaluator = Eval(
      jobMeta.workerPaths,
      Some(versionSupport.version),
      Vector.empty,
      jobMeta.fileResolver,
      Logger.Quiet
  )

  override val executorName: String = "dxExecutorWdl"

  override protected def typeAliases: Map[String, TSchema] = WdlUtils.toIRSchemaMap(wdlTypeAliases)

  override protected def evaluateInputs(
      jobInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    def getInputValues(inputTypes: Map[String, T]): Map[String, V] = {
      jobInputs.collect {
        case (name, (_, value)) if inputTypes.contains(name) =>
          name -> WdlUtils.fromIRValue(value, inputTypes(name), name)
        case (name, (t, v)) =>
          name -> WdlUtils.fromIRValue(v, WdlUtils.fromIRType(t), name)
      }
    }

    // This might be the input for the entire workflow or just a subblock.
    // If it is for a sublock, it may be for the body of a conditional or
    // scatter, in which case we only need the inputs of the body statements.
    val (inputTypes, inputValues) = jobMeta.blockPath match {
      case Vector() =>
        val inputTypes = workflow.inputs.map(inp => inp.name -> inp.wdlType).toMap
        // convert IR to WDL values
        val inputValues = getInputValues(inputTypes)
        if (logger.isVerbose) {
          logger.trace(s"""input parameters:
                          |${workflow.inputs
                            .map(TypeUtils.prettyFormatInput(_))
                            .mkString("\n")}
                          |input values:
                          |${WdlUtils.prettyFormatValues(inputValues)}""".stripMargin)
        }
        // evaluate - enable special handling for unset array values -
        // DNAnexus does not distinguish between null and empty for
        // array inputs, so we treat a null value for a non-optional
        // array that is allowed to be empty as the empty array.
        val inputBindings = InputOutput.inputsFromValues(workflow.name,
                                                         workflow.inputs,
                                                         inputValues,
                                                         evaluator,
                                                         ignoreDefaultEvalError = false,
                                                         nullCollectionAsEmpty = true)
        // convert back to IR
        (inputTypes, inputBindings.toMap)
      case path =>
        val block: WdlBlock =
          Block.getSubBlockAt(WdlBlock.createBlocks(workflow.body), path)
        val inputTypes = block.target match {
          case Some(conditional: TAT.Conditional) =>
            val (inputs, _) =
              WdlUtils.getClosureInputsAndOutputs(Vector(conditional), withField = true)
            inputs.map {
              case (name, (wdlType, _)) => name -> wdlType
            }
          case Some(scatter: TAT.Scatter) =>
            val (inputs, _) =
              WdlUtils.getClosureInputsAndOutputs(Vector(scatter), withField = true)
            inputs.map {
              case (name, (wdlType, _)) => name -> wdlType
            }
          case _ => block.inputs.map(inp => inp.name -> inp.wdlType).toMap
        }
        val inputValues = getInputValues(inputTypes)
        if (logger.isVerbose) {
          logger.trace(
              s"""input parameters:
                 |${inputTypes
                   .map {
                     case (name, wdlType) => s"  ${TypeUtils.prettyFormatType(wdlType)} ${name}"
                   }
                   .mkString("\n")}
                 |input values:
                 |${WdlUtils.prettyFormatValues(inputValues)}""".stripMargin
          )
        }
        (inputTypes, inputValues)
    }
    // convert back to IR
    inputValues.map {
      case (name, value) =>
        val wdlType = inputTypes(name)
        val irType = WdlUtils.toIRType(wdlType)
        val irValue = WdlUtils.toIRValue(value, wdlType)
        name -> (irType, irValue)
    }
  }

  override protected def evaluateOutputs(
      jobInputs: Map[String, (Type, Value)],
      addReorgStatus: Boolean
  ): Map[String, (Type, Value)] = {
    // This might be the output for the entire workflow or just a subblock.
    // If it is for a sublock, it may be for the body of a conditional or
    // scatter, in which case we only need the outputs of the body statements.
    val outputsParams = jobMeta.blockPath match {
      case Vector() => workflow.outputs
      case path =>
        val block: WdlBlock =
          Block.getSubBlockAt(WdlBlock.createBlocks(workflow.body), path)
        block.target match {
          case Some(conditional: TAT.Conditional) =>
            val (_, outputs) =
              WdlUtils.getClosureInputsAndOutputs(conditional.body, withField = true)
            outputs.values.toVector
          case Some(scatter: TAT.Scatter) =>
            val (_, outputs) =
              WdlUtils.getClosureInputsAndOutputs(scatter.body, withField = true)
            // exclude the scatter variable
            outputs.values.filter {
              case out: TAT.OutputParameter if out.name == scatter.identifier => false
              case _                                                          => true
            }.toVector
          case _ =>
            // we only need to expose the outputs that are not inputs
            val inputs = block.inputs.map(i => i.name -> i.wdlType).toMap
            block.outputs.filterNot(o => inputs.contains(o.name))
        }
    }
    val outputTypes = outputsParams.map(o => o.name -> o.wdlType).toMap
    // create the evaluation environment
    // first add the input values
    val inputWdlValues: Map[String, V] = jobInputs.map {
      case (name, (_, v)) if outputTypes.contains(name) =>
        name -> WdlUtils.fromIRValue(v, outputTypes(name), name)
      case (name, (t, v)) =>
        name -> WdlUtils.fromIRValue(v, WdlUtils.fromIRType(t), name)
    }
    // also add defaults for any optional values in the output
    // closure that are not inputs
    val outputWdlValues: Map[String, V] =
      WdlUtils.getOutputClosure(outputsParams).collect {
        case (name, T_Optional(_)) if !inputWdlValues.contains(name) =>
          // set missing optional inputs to null
          name -> V_Null
        case (name, T_Array(_, false)) if !inputWdlValues.contains(name) =>
          // DNAnexus does not distinguish between null and empty array inputs.
          // A non-optional, maybe-empty array may be missing, so set it to
          // the empty array.
          name -> V_Array()
      }
    if (logger.isVerbose) {
      logger.trace(
          s"""|input values: ${inputWdlValues}
              |output closure values: ${outputWdlValues} 
              |outputs parameters:
              |  ${outputsParams.map(TypeUtils.prettyFormatOutput(_)).mkString("\n  ")}
              |""".stripMargin
      )
    }
    // evaluate
    val env = WdlValueBindings(inputWdlValues ++ outputWdlValues)
    val evaluatedOutputValues = InputOutput.evaluateOutputs(outputsParams, evaluator, env)
    // convert back to IR
    val irOutputs = evaluatedOutputValues.toMap.map {
      case (name, value) =>
        val wdlType = outputTypes(name)
        val irType = WdlUtils.toIRType(wdlType)
        val irValue = WdlUtils.toIRValue(value, wdlType)
        name -> (irType, irValue)
    }
    if (addReorgStatus) {
      irOutputs + (Constants.ReorgStatus -> (TString, VString(Constants.ReorgStatusCompleted)))
    } else {
      irOutputs
    }
  }

  private def evaluateExpression(expr: TAT.Expr, wdlType: T, env: Map[String, (T, V)]): V = {
    evaluator.applyExprAndCoerce(expr, wdlType, Eval.createBindingsFromEnv(env))
  }

  private def getBlockOutputs(elements: Vector[TAT.WorkflowElement]): Map[String, T] = {
    val (_, outputs) = WdlUtils.getClosureInputsAndOutputs(elements, withField = false)
    outputs.values.map {
      case TAT.OutputParameter(name, wdlType, _, _) => name -> wdlType
    }.toMap
  }

  /**
    * Recursively evaluate PrivateVariables in WorkflowElements.
    * This method is exposed so that we can unit-test it.
    * @param elements WorkflowElements
    * @param env initial environment
    * @return values from all nested variables (not including the initial env)
    */
  private[wdl] def evaluateWorkflowElementVariables(
      elements: Seq[TAT.WorkflowElement],
      env: Map[String, (T, V)]
  ): Map[String, (T, V)] = {
    elements.foldLeft(Map.empty[String, (T, V)]) {
      case (accu, TAT.PrivateVariable(name, wdlType, expr, _)) =>
        val value = evaluateExpression(expr, wdlType, accu ++ env)
        accu + (name -> (wdlType, value))
      case (accu, TAT.Conditional(expr, body, _)) =>
        // evaluate the condition
        val results = evaluateExpression(expr, expr.wdlType, accu ++ env) match {
          case V_Boolean(true) =>
            // condition is true, evaluate the internal block.
            evaluateWorkflowElementVariables(body, accu ++ env).map {
              case (key, (t, value)) => key -> (T_Optional(t), V_Optional(value))
            }
          case V_Boolean(false) =>
            // condition is false, return V_Null for all the values
            getBlockOutputs(body).map {
              case (key, wdlType) => key -> (T_Optional(wdlType), V_Null)
            }
          case other =>
            throw new AppInternalException(s"Unexpected condition expression value ${other}")
        }
        accu ++ results
      case (accu, TAT.Scatter(id, expr, body, _)) =>
        val collection: Vector[V] =
          evaluateExpression(expr, expr.wdlType, accu ++ env) match {
            case V_Array(array) => array
            case other =>
              throw new AppInternalException(s"Unexpected class ${other.getClass}, ${other}")
          }
        val outputTypes: Map[String, T] = getBlockOutputs(body)
        val outputNames: Vector[String] = outputTypes.keys.toVector
        // iterate on the collection, evaluate the body N times,
        // transpose the results into M vectors of N items
        val outputValues: Vector[Vector[V]] =
          collection.map { v =>
            val envInner = accu ++ env + (id -> (expr.wdlType, v))
            val bodyValues = evaluateWorkflowElementVariables(body, envInner)
            outputNames.map(bodyValues(_)._2)
          }.transpose
        // Add the WDL array type to each vector
        val results = outputNames.zip(outputValues).map {
          case (name, values) =>
            val arrayType = T_Array(outputTypes(name), nonEmpty = false)
            val arrayValue = V_Array(values)
            name -> (arrayType, arrayValue)
        }
        accu ++ results
      case (_, other) =>
        throw new Exception(s"type ${other.getClass} while evaluating expressions")
    }
  }

  object WdlBlockContext {
    private[wdl] def evaluateCallInputs(
        call: TAT.Call,
        env: Map[String, (T, V)] = Map.empty
    ): Map[String, (T, V)] = {
      call.callee.input.flatMap {
        case (name, (wdlType, optional)) if call.inputs.contains(name) =>
          val optType = if (optional) {
            TypeUtils.ensureOptional(wdlType)
          } else {
            TypeUtils.unwrapOptional(wdlType)
          }
          val value = evaluateExpression(call.inputs(name), optType, env)
          Some(name -> (optType, value))
        case (name, (_, optional)) if optional =>
          logger.trace(s"no input for optional input ${name} to call ${call.fullyQualifiedName}")
          None
        case (name, _) =>
          throw new Exception(
              s"missing non-optional input ${name} to call ${call.fullyQualifiedName}"
          )
      }
    }
  }

  case class WdlBlockContext(block: WdlBlock, env: Map[String, (T, V)])
      extends BlockContext[WdlBlock] {
    private def call: TAT.Call = block.call
    private def dxApi = jobMeta.dxApi

    private def evaluateCallInputs(
        extraEnv: Map[String, (T, V)] = Map.empty
    ): Map[String, (T, V)] = {
      WdlBlockContext.evaluateCallInputs(call, env ++ extraEnv)
    }

    private def launchCall(
        callInputs: Map[String, (T, V)],
        nameDetail: Option[String] = None,
        folder: Option[String] = None
    ): (DxExecution, ExecutableLink, String) = {
      logger.traceLimited(
          s"""|call = ${call}
              |callInputs = ${callInputs}
              |""".stripMargin,
          minLevel = TraceLevel.VVerbose
      )
      val executableLink = getExecutableLink(call.callee.name)
      val callInputsIR = WdlUtils.toIR(callInputs)
      val instanceType = tasks.get(call.callee.name).flatMap { task =>
        val meta: Meta = Meta.create(versionSupport.version, task.meta)
        val isNative = meta.get("type", Vector(T_Boolean)) match {
          case Some(V_Boolean(b)) => b
          case _                  => false
        }
        if (isNative) {
          None
        } else {
          val callIO = TaskInputOutput(task, logger)
          val inputWdlValues: Map[String, V] = callInputsIR.collect {
            case (name, (t, v)) if !name.endsWith(ParameterLink.FlatFilesSuffix) =>
              val wdlType = WdlUtils.fromIRType(t, wdlTypeAliases)
              name -> WdlUtils.fromIRValue(v, wdlType, name)
          }
          // add default values for any missing inputs
          val callInputs =
            callIO.inputsFromValues(inputWdlValues, evaluator, ignoreDefaultEvalError = false)
          val runtime =
            Runtime(versionSupport.version,
                    task.runtime,
                    task.hints,
                    evaluator,
                    None,
                    ctx = Some(callInputs))
          try {
            val request = runtime.parseInstanceType
            val instanceType = jobMeta.instanceTypeDb.apply(request)
            logger.traceLimited(s"Precalculated instance type for ${task.name}: ${instanceType}")
            Some(instanceType)
          } catch {
            case e: Throwable =>
              logger.traceLimited(
                  s"""|Failed to precalculate the instance type for
                      |task ${task.name}.
                      |
                      |${e}
                      |""".stripMargin
              )
              None
          }
        }
      }

      val (dxExecution, execName) =
        launchJob(executableLink,
                  call.actualName,
                  callInputsIR,
                  nameDetail,
                  instanceType.map(_.name),
                  folder = folder,
                  prefixOutputs = true)
      (dxExecution, executableLink, execName)
    }

    private def launchCall(blockIndex: Int): Map[String, ParameterLink] = {
      val callInputs = evaluateCallInputs()
      val (dxExecution, executableLink, callName) =
        launchCall(callInputs, folder = Some(blockIndex.toString))
      jobMeta.createExecutionOutputLinks(dxExecution, executableLink.outputs, Some(callName))
    }

    private val qualifiedNameRegexp = "(.+)\\.(.+)".r

    private def lookup(name: String, env: Map[String, (T, V)]): Option[V] = {
      def inner(innerName: String): Option[V] = {
        innerName match {
          case _ if env.contains(innerName) =>
            Some(env(innerName)._2)
          case qualifiedNameRegexp(lhs, rhs) =>
            inner(lhs).map {
              case V_Pair(left, _) if rhs == "left" =>
                left
              case V_Pair(_, right) if rhs == "right" =>
                right
              case V_Struct(_, members) if members contains rhs =>
                members(rhs)
              case V_Call(_, members) if members contains rhs =>
                members(rhs)
              case V_Object(members) if members contains rhs =>
                members(rhs)
              case _ =>
                throw new Exception(s"field ${rhs} does not exist in ${lhs}")
            }
          case _ =>
            None
        }
      }
      inner(name)
    }

    private def prepareSubworkflowInputs(
        executableLink: ExecutableLink,
        extraEnv: Map[String, (T, V)] = Map.empty
    ): Map[String, (Type, Value)] = {
      val inputEnv = env ++ extraEnv
      logger.traceLimited(
          s"""|prepareSubworkflowInputs (${executableLink.name})
              |env:
              |${inputEnv.mkString("\n")}
              |
              |linkInfo: ${executableLink}
              |""".stripMargin,
          minLevel = TraceLevel.VVerbose
      )
      // Missing inputs may be optional or have a default value. If they are
      // actually missing, it will result in a platform error.
      executableLink.inputs.flatMap {
        case (name, t) =>
          val value = lookup(name, inputEnv)
          logger.traceLimited(s"lookupInEnv(${name} = ${value})", minLevel = TraceLevel.VVerbose)
          value.map { value =>
            val wdlType = WdlUtils.fromIRType(t, wdlTypeAliases)
            name -> (t, WdlUtils.toIRValue(value, wdlType))
          }
      }
    }

    private def launchConditional(): Map[String, ParameterLink] = {
      val cond = block.target match {
        case Some(TAT.Conditional(expr, _, _)) =>
          evaluateExpression(expr, T_Boolean, env)
        case _ =>
          throw new Exception(s"invalid conditional block ${block}")
      }
      val links = (cond, block.kind) match {
        case (V_Boolean(true), BlockKind.ConditionalOneCall) =>
          // A subblock containing exactly one call. For example:
          // if (flag) {
          //     call zinc as inc3 { input: a = num}
          // }
          // The flag evaluates to true, so execute the inner call
          // and ensure it's output type is optional
          launchCall(block.index)
        case (V_Boolean(true), BlockKind.ConditionalComplex) =>
          // A complex subblock requiring a fragment runner, or a subworkflow. For example:
          //
          // if (flag) {
          //    call zinc as inc3 { input: a = num}
          //    call zinc as inc4 { input: a = num + 3 }
          //
          //    Int b = inc4.result + 14
          //    call zinc as inc5 { input: a = b * 4 }
          // }
          //
          // There must be exactly one sub-workflow.
          assert(execLinkInfo.size == 1)
          val executableLink = execLinkInfo.values.head
          val callInputs = prepareSubworkflowInputs(executableLink)
          // there is no good human-readable name for a conditional, so we just
          // use the current block index
          val (dxExecution, _) =
            launchJob(executableLink,
                      executableLink.name,
                      callInputs,
                      folder = Some(block.index.toString))
          jobMeta.createExecutionOutputLinks(dxExecution, executableLink.outputs)
        case (V_Boolean(false), _) => Map.empty
        case _ =>
          throw new Exception(s"invalid conditional value ${cond}")
      }
      // Add optional modifier to the return types.
      links.map {
        case (key, link) => key -> link.makeOptional
      }
    }

    /*
    Scatter jobs are implemented as independent jobs for each element of
    the scatter, and a collection job that depends on the completion of
    all the scatter element jobs. This is necessary when the output is a
    non-native DNAx type. For example, the math workflow below calls a
    scatter where each job returns an array of files. The GenFiles.result
    is a ragged array of files (Array[Array[File]]). The conversion between
    these two types is difficult, and requires this applet.
    ```
    task GenFiles {
      ...
      output {
          Array[File] result
      }
    }
    workflow math {
        scatter (k in [2,3,5]) {
            call GenFiles { input: len=k }
        }
        output {
            Array[Array[File] result = GenFiles.result
        }
    }
    ```
    Diagram
              scatter
             /   | .. \
       child-jobs      \
                        \
                         collect
    Design
      The collect applet takes three inputs:
    1) job-ids     (array of strings)
    2) field names (array of strings)
    3) WDL types   (array of strings)
      It waits for all the scatter child jobs to complete, using the dependsOn field.
    For each field F:
      - Get the value of F from all the child jobs
      - Merge. This is complex but doable for non native dx types. For
        example, to merge the GenFiles output, we need to merge an array
        of array of files into a hash with a companion flat array of
        files.
    outputs: the merged value for each field.
    Larger context
      The parent scatter returns ebors to each of the collect output fields. This
    allows it to return immediately, and not wait for the child jobs to complete.
    Each scatter requires its own collect applet, because the output type is
    the same as the scatter output type.

    Continuations
      A scatter may contain many sub-jobs. Rather than enforce a maximum number of
      scatter sub-jobs, we instead chain scatters so that any number of sub-jobs
      can be executed with a maximum number running at one time. For example:

      scatter (i in range(2000)) { ... }

      translates to:

      scatter(1-1000, jobId=A, parents=[])
      |_exec job 1..1000
      |_exec scatter(1001..2000, jobId=B, parents=[A]) // does not run until jobs 1-1000 are complete
             |_exec job 1001..2000
             |_exec collect(parents=[A,B]) // does not run until jobs 1001-2000 are complete
     */

    private def evaluateScatterCollection(expr: TAT.Expr): (Vector[V], Option[Int]) = {
      evaluateExpression(expr, expr.wdlType, env) match {
        case V_Array(array) if jobMeta.scatterStart == 0 && array.size <= jobMeta.scatterSize =>
          (array, None)
        case V_Array(array) =>
          val scatterEnd = jobMeta.scatterStart + jobMeta.scatterSize
          if (scatterEnd < array.size) {
            (array.slice(jobMeta.scatterStart, scatterEnd), Some(scatterEnd))
          } else {
            (array.drop(jobMeta.scatterStart), None)
          }
        case other =>
          throw new AppInternalException(s"scatter value ${other} is not an array")
      }
    }

    private def truncate(scatterName: String): String = {
      if (scatterName.length > WorkflowExecutor.JobNameLengthLimit) {
        s"${scatterName.substring(0, WorkflowExecutor.JobNameLengthLimit - 3)}..."
      } else {
        scatterName
      }
    }

    private def getFileName(path: String): String = {
      jobMeta.fileResolver.resolve(path) match {
        case local: LocalFileSource => local.canonicalPath.getFileName.toString
        case dx: DxFileSource       => dx.dxFile.getName
        case other                  => other.toString
      }
    }

    // create a short, easy to read, description for a scatter element.
    private[wdl] def getScatterName(item: V): Option[String] = {
      item match {
        case _ if EvalUtils.isPrimitive(item) => Some(truncate(EvalUtils.formatPrimitive(item)))
        case V_File(path)                     => Some(truncate(getFileName(path)))
        case V_Directory(path)                => Some(truncate(getFileName(path)))
        case V_Optional(x)                    => getScatterName(x)
        case V_Pair(l, r) =>
          val ls = getScatterName(l)
          val rs = getScatterName(r)
          (ls, rs) match {
            case (Some(ls1), Some(rs1)) => Some(truncate(s"(${ls1}, ${rs1})"))
            case _                      => None
          }
        case V_Array(array) =>
          val itemStr =
            WdlWorkflowExecutor.getComplexScatterName(array.iterator.map(getScatterName))
          Some(s"[${itemStr}]")
        case V_Map(members) =>
          val memberStr = WdlWorkflowExecutor.getComplexScatterName(
              members.iterator.map {
                case (k, v) =>
                  (getScatterName(k), getScatterName(v)) match {
                    case (Some(keyStr), Some(valStr)) => Some(s"${keyStr}: ${valStr}")
                    case _                            => None
                  }
              }
          )
          Some(s"{${memberStr}}")
        case V_Object(members) =>
          val memberStr = WdlWorkflowExecutor.getComplexScatterName(
              members.iterator.map {
                case (k, v) =>
                  getScatterName(v) match {
                    case Some(valStr) => Some(s"${k}: ${valStr}")
                    case _            => None
                  }
              }
          )
          Some(s"{${memberStr}}")
        case V_Struct(name, members) =>
          val memberStr = WdlWorkflowExecutor.getComplexScatterName(
              members.iterator.map {
                case (k, v) =>
                  getScatterName(v) match {
                    case Some(valStr) => Some(s"${k}: ${valStr}")
                    case _            => None
                  }
              }
          )
          Some(s"${name} ${memberStr}")
        case _ =>
          None
      }
    }

    private def launchScatterCallJobs(identifier: String,
                                      itemType: T,
                                      collection: Vector[V]): Vector[DxExecution] = {
      collection.map { item =>
        val callInputs = evaluateCallInputs(Map(identifier -> (itemType, item)))
        val callNameDetail = getScatterName(item)
        val (dxExecution, _, _) = launchCall(callInputs, callNameDetail)
        dxExecution
      }
    }

    private def launchScatterSubblockJobs(identifier: String,
                                          itemType: T,
                                          collection: Vector[V]): Vector[DxExecution] = {
      assert(execLinkInfo.size == 1)
      val executableLink = execLinkInfo.values.head
      collection.map { item =>
        val callInputs =
          prepareSubworkflowInputs(executableLink, Map(identifier -> (itemType, item)))
        val callNameDetail = getScatterName(item)
        val (dxExecution, _) =
          launchJob(executableLink, executableLink.name, callInputs, callNameDetail)
        dxExecution
      }
    }

    private def prepareBlockOutputs(
        outputs: Map[String, ParameterLink]
    ): Map[String, ParameterLink] = {
      if (jobMeta.useManifests && outputs.contains(Constants.OutputManifest)) {
        outputs
      } else {
        val outputNames = block.outputNames
        logger.traceLimited(
            s"""|processOutputs
                |  env = ${env.keys}
                |  fragResults = ${outputs.keys}
                |  exportedVars = ${outputNames}
                |""".stripMargin,
            minLevel = TraceLevel.VVerbose
        )
        val inputsIR = WdlUtils.toIR(env.view.filterKeys(outputNames.contains).toMap)
        val inputLinks = jobMeta.createOutputLinks(inputsIR, validate = false)
        val outputLink = outputs.view.filterKeys(outputNames.contains).toMap
        inputLinks ++ outputLink
      }
    }

    /**
      * stick the IDs of all the parent jobs and all the child jobs to exclude
      * (i.e. the continue/collect jobs) into details - we'll use these in the
      * collect step
      * @param nextStart index at which to start the next scatter
      * @return
      */
    private def createSubjobDetails(nextStart: Option[Int] = None): JsValue = {
      val parents = jobMeta.getJobDetail(WorkflowExecutor.ParentsKey) match {
        case Some(JsArray(array)) => array.map(JsUtils.getString(_))
        case _                    => Vector.empty
      }
      // add the current job to the list of parents
      val allParents = parents :+ jobMeta.jobId
      val details = Map(WorkflowExecutor.ParentsKey -> JsArray(allParents.map(JsString(_))))
      val continueDetails = nextStart match {
        case Some(i) => Map(Constants.ContinueStart -> JsNumber(i))
        case None    => Map.empty
      }
      JsObject(details ++ continueDetails)
    }

    private def prepareScatterResults(dxSubJob: DxExecution): Map[String, ParameterLink] = {
      val resultTypes: Map[String, Type] = block.outputs.map {
        case TAT.OutputParameter(name, wdlType, _, _) =>
          name -> WdlUtils.toIRType(wdlType)
      }.toMap
      // Return JBORs for all the outputs. Since the signature of the sub-job
      // is exactly the same as the parent, we can immediately exit the parent job.
      val links = jobMeta.createExecutionOutputLinks(dxSubJob, resultTypes)
      if (logger.isVerbose) {
        val linkStr = links.mkString("\n")
        logger.traceLimited(s"resultTypes=${resultTypes}")
        logger.traceLimited(s"promises=${linkStr}")
      }
      links
    }

    /**
      * Lauch a job to continue a large scatter.
      * @param childJobs child jobs on which the continue job will depend
      * @param nextStart the index at which to continue the scatter
      * @return
      */
    private def launchScatterContinue(
        childJobs: Vector[DxExecution],
        nextStart: Int
    ): Map[String, ParameterLink] = {
      assert(childJobs.nonEmpty)
      // Run a sub-job with the "continue" entry point.
      // We need to provide the exact same inputs.
      val dxSubJob: DxExecution = dxApi.runSubJob(
          "continue",
          Some(jobMeta.instanceTypeDb.defaultInstanceType.name),
          JsObject(jobMeta.rawJsInputs),
          childJobs,
          jobMeta.delayWorkspaceDestruction,
          Some(s"continue_scatter($nextStart)"),
          Some(createSubjobDetails(Some(nextStart)))
      )
      prepareScatterResults(dxSubJob)
    }

    private def launchScatterCollect(childJobs: Vector[DxExecution]): Map[String, ParameterLink] = {
      assert(childJobs.nonEmpty)
      // Run a sub-job with the "collect" entry point.
      // We need to provide the exact same inputs.
      val dxSubJob: DxExecution = dxApi.runSubJob(
          "collect",
          Some(jobMeta.instanceTypeDb.defaultInstanceType.name),
          JsObject(jobMeta.rawJsInputs),
          childJobs,
          jobMeta.delayWorkspaceDestruction,
          Some(s"collect_scatter"),
          Some(createSubjobDetails())
      )
      prepareScatterResults(dxSubJob)
    }

    private def launchScatter(): Map[String, ParameterLink] = {
      val (identifier, itemType, collection, next) = block.target match {
        case Some(TAT.Scatter(identifier, expr, _, _)) =>
          val (collection, next) = evaluateScatterCollection(expr)
          val itemType = expr.wdlType match {
            case T_Array(t, _) => t
            case _ =>
              throw new Exception(s"scatter type ${expr.wdlType} is not an array")
          }
          (identifier, itemType, collection, next)
        case _ =>
          throw new RuntimeException(s"invalid scatter block ${block}")
      }
      val childJobs: Vector[DxExecution] = block.kind match {
        case BlockKind.ScatterOneCall =>
          launchScatterCallJobs(identifier, itemType, collection)
        case BlockKind.ScatterComplex =>
          launchScatterSubblockJobs(identifier, itemType, collection)
        case _ =>
          throw new RuntimeException(s"invalid scatter block ${block}")
      }
      next match {
        case Some(index) =>
          // there are remaining chunks - call a continue sub-job
          launchScatterContinue(childJobs, index)
        case None =>
          // this is the last chunk - call collect sub-job to gather all the results
          launchScatterCollect(childJobs)
      }
    }

    private case class ChildExecution(execName: String,
                                      seqNum: Int,
                                      outputs: Map[String, JsValue],
                                      exec: DxExecution)

    private def parseOneResult(value: JsValue, excludeIds: Set[String]): Option[ChildExecution] = {
      val fields = value.asJsObject.fields
      val (exec, desc) = fields.get("id") match {
        case Some(JsString(id)) if excludeIds.contains(id) =>
          logger.trace(s"Ignoring result for job ${id}")
          return None
        case Some(JsString(id)) if id.startsWith("job-") =>
          val job = dxApi.job(id)
          val desc = fields("describe").asJsObject
          (job, desc)
        case Some(JsString(id)) if id.startsWith("analysis-") =>
          val analysis = dxApi.analysis(id)
          val desc = fields("describe").asJsObject
          (analysis, desc)
        case Some(other) =>
          throw new Exception(s"malformed id field ${other.prettyPrint}")
        case None =>
          throw new Exception(s"field id not found in ${value.prettyPrint}")
      }
      logger.trace(s"parsing desc ${desc} for ${exec}")
      val (execName, details, output) =
        desc.getFields("executableName", "details", "output") match {
          case Seq(JsString(execName), JsObject(details), JsObject(output)) =>
            (execName, details, output)
        }
      val seqNum = details.get(WorkflowExecutor.SeqNumber) match {
        case Some(JsNumber(i)) => i.toIntExact
        case other             => throw new Exception(s"Invalid seqNumber ${other}")
      }
      Some(ChildExecution(execName, seqNum, output, exec))
    }

    private def submitRequest(
        parentJobId: Option[String],
        cursor: JsValue,
        excludeIds: Set[String],
        limit: Option[Int]
    ): (Vector[ChildExecution], JsValue) = {
      val parentField: Map[String, JsValue] = parentJobId match {
        case None     => Map.empty
        case Some(id) => Map("parentJob" -> JsString(id))
      }
      val cursorField: Map[String, JsValue] = cursor match {
        case JsNull      => Map.empty
        case cursorValue => Map("starting" -> cursorValue)
      }
      val limitField: Map[String, JsValue] = limit match {
        case None    => Map.empty
        case Some(i) => Map("limit" -> JsNumber(i))
      }
      val describeField: Map[String, JsValue] = Map(
          "describe" -> JsObject(
              "fields" -> DxObject
                .requestFields(Set(Field.Output, Field.ExecutableName, Field.Details))
          )
      )
      val response = dxApi.findExecutions(parentField ++ cursorField ++ limitField ++ describeField)
      val results: Vector[ChildExecution] =
        response.fields.get("results") match {
          case Some(JsArray(results)) =>
            results.flatMap(res => parseOneResult(res, excludeIds))
          case Some(other) =>
            throw new Exception(s"malformed results field ${other.prettyPrint}")
          case None =>
            throw new Exception(s"missing results field ${response}")
        }
      (results, response.fields("next"))
    }

    private def findChildExecutions(parentJobId: Option[String],
                                    excludeIds: Set[String],
                                    limit: Option[Int] = None): Vector[ChildExecution] = {
      Iterator
        .unfold[Vector[ChildExecution], Option[JsValue]](Some(JsNull)) {
          case None => None
          case Some(cursor: JsValue) =>
            submitRequest(parentJobId, cursor, excludeIds, limit) match {
              case (Vector(), _)     => None
              case (results, JsNull) => Some(results, None)
              case (results, next)   => Some(results, Some(next))
            }
        }
        .toVector
        .flatten
        .sortWith(_.seqNum < _.seqNum)
    }

    /**
      * Gets all the jobs launched by this job's origin job, excluding
      * any continue and collect sub-jobs.
      * @return
      */
    private def getScatterJobs: Vector[ChildExecution] = {
      val childExecs = jobMeta.getJobDetail(WorkflowExecutor.ParentsKey) match {
        case Some(JsArray(array)) =>
          val parentJobIds = array.map(JsUtils.getString(_))
          val excludeJobIds = parentJobIds.toSet + jobMeta.jobId
          parentJobIds.flatMap { parentJobId =>
            findChildExecutions(Some(parentJobId), excludeJobIds)
          }
        case _ =>
          val parentJob = jobMeta.parentJob match {
            case Some(job) => job
            case None =>
              throw new Exception(s"Can't get parent job for $jobMeta.jobDesc")
          }
          findChildExecutions(Some(parentJob.id), Set(jobMeta.jobId))
      }
      logger.trace(s"childExecs=${childExecs}")
      childExecs
    }

    private def collectScatter(): Map[String, ParameterLink] = {
      val childExecutions = getScatterJobs

      val outputTypes: Map[String, (String, Type)] = block.kind match {
        case BlockKind.ScatterOneCall =>
          call.callee.output.map {
            case (name, wdlType) =>
              val fqn = s"${call.actualName}.${name}"
              val irType = WdlUtils.toIRType(wdlType)
              fqn -> (name, irType)
          }
        case BlockKind.ScatterComplex =>
          assert(execLinkInfo.size == 1)
          execLinkInfo.values.head.outputs.map {
            case (name, irType) => name -> (name, irType)
          }
        case _ =>
          throw new RuntimeException(s"invalid block ${block}")
      }

      val (manifestId, executionName, childOutputs) =
        if (jobMeta.useManifests) {
          // each job has an output manifest - we need to download them all
          childExecutions
            .foldLeft(Option.empty[String],
                      Option.empty[String],
                      Vector.empty[Map[String, JsValue]]) {
              case ((id, execName, childOutputs), childExec) =>
                val manifestFile = childExec.outputs.get(Constants.OutputManifest) match {
                  case Some(fileObj: JsObject) if DxFile.isLinkJson(fileObj) =>
                    Some(DxFile.fromJson(dxApi, fileObj))
                  case Some(JsString(uri)) if uri.startsWith(DxPath.DxUriPrefix) =>
                    Some(dxApi.resolveFile(uri))
                  case None =>
                    // maybe the applet doesn't have any outputs
                    None
                  case other =>
                    throw new Exception(s"invalid manifest file value ${other}")
                }
                val (manifestId, manifestValues) = manifestFile
                  .map { dxFile =>
                    val manifestJson = new String(dxApi.downloadBytes(dxFile)).parseJson
                    val manifest = Manifest.parse(manifestJson)
                    (manifest.id, manifest.jsValues)
                  }
                  .getOrElse((None, Map.empty[String, JsValue]))
                // all scatter jobs manifests should have the same ID
                val newId = (id, manifestId) match {
                  case (None, None)                         => None
                  case (None, Some(id))                     => Some(id)
                  case (Some(id), None)                     => Some(id)
                  case (Some(id1), Some(id2)) if id1 == id2 => Some(id1)
                  case (Some(id1), Some(id2)) =>
                    throw new Exception(
                        s"scatter job output manifests had different IDs: ${id1} != ${id2}"
                    )
                }
                val newExecName = (execName, childExec.execName) match {
                  case (None, execName)                                       => Some(execName)
                  case (Some(execName1), execName2) if execName1 == execName2 => Some(execName1)
                  case (Some(execName1), execName2) =>
                    throw new Exception(
                        s"child execution names do not match: ${execName1} != ${execName2}"
                    )
                }
                (newId, newExecName, childOutputs :+ manifestValues)
            }
        } else {
          (None, None, childExecutions.map(_.outputs))
        }

      val arrayValues: Map[String, (Type, Value)] = outputTypes.view.mapValues {
        case (name, irType) =>
          val arrayType = TArray(irType)
          val nameEncoded = Parameter.encodeDots(name)
          val longNameEncoded = executionName.map(e => Parameter.encodeDots(s"${e}.${name}"))
          val arrayValue = childOutputs.flatMap { outputs =>
            val jsValue = outputs.get(nameEncoded).orElse(longNameEncoded.flatMap(outputs.get))
            (irType, jsValue) match {
              case (_, Some(jsValue)) =>
                Some(jobMeta.inputDeserializer.deserializeInputWithType(jsValue, irType))
              case (TOptional(_), None) => None
              case (_, None)            =>
                // Required output that is missing
                throw new Exception(s"missing required field <${name}> in results")
            }
          }
          (arrayType, VArray(arrayValue))
      }.toMap

      if (arrayValues.isEmpty) {
        Map.empty
      } else if (jobMeta.useManifests) {
        if (manifestId.isEmpty) {
          throw new Exception("missing manifest Id")
        }
        // upload the merged manifest file
        val outputJson = arrayValues.map {
          case (name, (t, v)) => Parameter.encodeDots(name) -> ValueSerde.serializeWithType(v, t)
        }
        val manifest = Manifest(outputJson, id = manifestId)
        val destination = s"${jobMeta.manifestFolder}/${jobMeta.jobId}_output.manifest.json"
        val manifestDxFile = dxApi.uploadString(manifest.toJson.prettyPrint, destination)
        val outputValues = Map(
            Constants.OutputManifest -> (TFile, VFile(manifestDxFile.asUri))
        )
        jobMeta.createOutputLinks(outputValues, validate = false)
      } else {
        jobMeta.createOutputLinks(arrayValues, validate = false)
      }
    }

    override def launch(): Map[String, ParameterLink] = {
      val outputs: Map[String, ParameterLink] = block.kind match {
        case BlockKind.CallWithSubexpressions | BlockKind.CallFragment =>
          launchCall(block.index)
        case BlockKind.ConditionalOneCall | BlockKind.ConditionalComplex =>
          launchConditional()
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          assert(jobMeta.scatterStart == 0)
          launchScatter()
        case BlockKind.ExpressionsOnly =>
          Map.empty
        case BlockKind.CallDirect =>
          throw new RuntimeException("unreachable state")
      }
      prepareBlockOutputs(outputs)
    }

    override def continue(): Map[String, ParameterLink] = {
      val outputs: Map[String, ParameterLink] = block.kind match {
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          assert(jobMeta.scatterStart > 0)
          launchScatter()
        case _ =>
          throw new RuntimeException(s"cannot continue non-scatter block ${block}")
      }
      prepareBlockOutputs(outputs)
    }

    override def collect(): Map[String, ParameterLink] = {
      val outputs: Map[String, ParameterLink] = block.kind match {
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          collectScatter()
        case _ =>
          throw new RuntimeException(s"cannot continue non-scatter block ${block}")
      }
      prepareBlockOutputs(outputs)
    }

    override lazy val prettyFormat: String = {
      val envStr = if (env.isEmpty) {
        " <empty>"
      } else {
        s"\n${WdlUtils.prettyFormatEnv(env)}"
      }
      s"""${block.prettyFormat}
         |Env:${envStr}
         |""".stripMargin
    }
  }

  private val compoundNameRegexp = ".+[.\\[].+".r

  override def evaluateBlockInputs(
      jobInputs: Map[String, (Type, Value)]
  ): WdlBlockContext = {
    val block: WdlBlock =
      Block.getSubBlockAt(WdlBlock.createBlocks(workflow.body), jobMeta.blockPath)
    // convert the inputs to WDL
    val initEnv = jobInputs.map {
      case (name, (irType, irValue)) =>
        val wdlType = WdlUtils.fromIRType(irType, wdlTypeAliases)
        val wdlValue = WdlUtils.fromIRValue(irValue, wdlType, name)
        name -> (wdlType, wdlValue)
    }
    val inputEnv = block.inputs.foldLeft(initEnv) {
      case (accu, blockInput: WdlBlockInput) if accu.contains(blockInput.name) =>
        accu
      case (accu, _: ComputedBlockInput) =>
        // this is the scatter variable, or some expression that references the scatter
        // variable - we can ignore it since it will be evaluated during launching of
        // the scatter jobs
        accu
      case (accu, RequiredBlockInput(name, wdlType)) if compoundNameRegexp.matches(name) =>
        // the input name is a compound reference - evaluate it as an identifier
        val expr = versionSupport.parseExpression(name, DefaultBindings(accu.map {
          case (name, (wdlType, _)) => name -> wdlType
        }), docSource)
        accu + (name -> (wdlType, evaluateExpression(expr, wdlType, accu)))
      case (_, RequiredBlockInput(name, _)) =>
        throw new Exception(s"missing required input ${name}")
      case (accu, OverridableBlockInputWithStaticDefault(name, wdlType, defaultValue)) =>
        accu + (name -> (wdlType, defaultValue))
      case (accu, OverridableBlockInputWithDynamicDefault(name, wdlType, defaultExpr)) =>
        val wdlValue = evaluateExpression(defaultExpr, wdlType, accu)
        accu + (name -> (wdlType, wdlValue))
      case (accu, OptionalBlockInput(name, wdlType)) =>
        // the input is missing but it could be optional - add a V_Null
        accu + (name -> (wdlType, V_Null))
    }
    val prereqEnv = evaluateWorkflowElementVariables(block.prerequisites, inputEnv)
    WdlBlockContext(block, inputEnv ++ prereqEnv)
  }
}
