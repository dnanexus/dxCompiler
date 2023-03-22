package dx.executor.wdl

import dx.AppInternalException
import dx.api.DxExecution
import dx.core.Constants
import dx.core.ir.{Block, BlockKind, DxName, ExecutableLink, ParameterLink, Type, Value}
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
  WdlBlockOutput,
  WdlDxName,
  WdlOptions,
  WdlUtils
}
import dx.executor.{JobMeta, WorkflowExecutor}
import dx.util.{DefaultBindings, FileNode, TraceLevel}
import dx.util.exceptionToString
import spray.json.JsValue
import wdlTools.eval.{Eval, EvalException, EvalUtils, WdlValueBindings}
import wdlTools.eval.WdlValues._
import wdlTools.exec.{InputOutput, TaskInputOutput}
import wdlTools.types.{TypeUtils, TypedAbstractSyntax => TAT}
import wdlTools.types.WdlTypes._

object WdlWorkflowExecutor {
  def create(jobMeta: JobMeta, separateOutputs: Boolean): WdlWorkflowExecutor = {
    // parse the workflow source code to get the WDL document
    val wdlOptions = jobMeta.parserOptions.map(WdlOptions.fromJson).getOrElse(WdlOptions.default)
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceString(jobMeta.sourceCode, wdlOptions, jobMeta.fileResolver)
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
}

// TODO: implement graph building from workflow - input and output parameters should be ordered
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
    extends WorkflowExecutor[WdlBlock](jobMeta = jobMeta, separateOutputs = separateOutputs) {
  private val logger = jobMeta.logger
  private lazy val evaluator = Eval(
      jobMeta.workerPaths,
      Some(versionSupport.version),
      Vector.empty,
      jobMeta.fileResolver
  )

  override val executorName: String = "dxExecutorWdl"

  override protected lazy val typeAliases: Map[String, TSchema] =
    WdlUtils.toIRSchemaMap(wdlTypeAliases)

  override protected def evaluateInputs(
      jobInputs: Map[DxName, (Type, Value)]
  ): Map[DxName, (Type, Value)] = {
    def getInputValues(inputTypes: Map[DxName, T]): Map[DxName, V] = {
      jobInputs.collect {
        case (name, (_, value)) if inputTypes.contains(name) =>
          name -> WdlUtils.fromIRValue(value, inputTypes(name), name.decoded)
        case (name, (t, v)) =>
          name -> WdlUtils.fromIRValue(v, WdlUtils.fromIRType(t, wdlTypeAliases), name.decoded)
      }
    }

    // This might be the input for the entire workflow or just a subblock.
    // If it is for a sublock, it may be for the body of a conditional or
    // scatter, in which case we only need the inputs of the body statements.
    val (inputTypes: Map[DxName, T], inputValues: Map[DxName, V]) =
      jobMeta.blockPath match {
        case Vector() =>
          val inputTypes: Map[DxName, T] = workflow.inputs
            .map(inp => WdlDxName.fromSourceName(inp.name) -> inp.wdlType)
            .toMap
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
          val inputBindings =
            InputOutput.inputsFromValues(workflow.name, workflow.inputs, inputValues.map {
              case (dxName, v) => dxName.decoded -> v
            }, evaluator, ignoreDefaultEvalError = false, nullCollectionAsEmpty = true)
          (inputTypes, inputBindings.toMap.map {
            case (name, v) => WdlDxName.fromSourceName(name) -> v
          })
        case path =>
          val block: WdlBlock =
            Block.getSubBlockAt(WdlBlock.createBlocks(workflow.body), path)
          val targetElement = block.target match {
            case Some(conditional: TAT.Conditional) => Some(conditional)
            case Some(scatter: TAT.Scatter)         => Some(scatter)
            case _                                  => None
          }
          val targetInputTypes = targetElement
            .map { e =>
              val (inputs, _) =
                WdlUtils.getClosureInputsAndOutputs(Vector(e), withField = true)
              inputs.map {
                case (dxName, (wdlType, _)) => dxName -> wdlType
              }.toMap
            }
            .getOrElse(Map.empty[DxName, T])
          val exprInputTypes = block.inputs
            .filterNot(i => targetInputTypes.contains(i.name))
            .map(inp => inp.name -> inp.wdlType)
            .toMap
          val inputTypes = targetInputTypes ++ exprInputTypes
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
      jobInputs: Map[DxName, (Type, Value)],
      addReorgStatus: Boolean
  ): Map[DxName, (Type, Value)] = {
    // This might be the output for the entire workflow or just a subblock.
    // If it is for a sublock, it may be for the body of a conditional or
    // scatter, in which case we only need the outputs of the body statements.
    val outputParams: Vector[(DxName, (T, TAT.Expr))] = jobMeta.blockPath match {
      case Vector() =>
        workflow.outputs.map {
          case TAT.OutputParameter(name, wdlType, expr) =>
            WdlDxName.fromSourceName(name) -> (wdlType, expr)
        }
      case path =>
        val block: WdlBlock =
          Block.getSubBlockAt(WdlBlock.createBlocks(workflow.body), path)
        block.target match {
          case Some(conditional: TAT.Conditional) =>
            val (_, outputs) =
              WdlUtils.getClosureInputsAndOutputs(conditional.body, withField = true)
            outputs
          case Some(scatter: TAT.Scatter) =>
            val (_, outputs) =
              WdlUtils.getClosureInputsAndOutputs(scatter.body, withField = true)
            // exclude the scatter variable
            outputs.filterNot(_._1.decoded == scatter.identifier)
          case _ =>
            // we only need to expose the outputs that are not inputs
            val inputNames = block.inputs.map(_.name).toSet
            block.outputs
              .filterNot(o => inputNames.contains(o.name))
              .map {
                case WdlBlockOutput(name, wdlType, expr) => name -> (wdlType, expr)
              }
        }
    }
    val outputParamMap = outputParams.toMap
    // create the evaluation environment
    // first add the input values
    val inputWdlValues: Map[DxName, V] = jobInputs.map {
      case (dxName, (_, v)) if outputParamMap.contains(dxName) =>
        dxName -> WdlUtils.fromIRValue(v, outputParamMap(dxName)._1, dxName.decoded)
      case (dxName, (t, v)) =>
        dxName -> WdlUtils.fromIRValue(v, WdlUtils.fromIRType(t, wdlTypeAliases), dxName.decoded)
    }
    // next add defaults for any optional values in the output closure that are not inputs
    val outputWdlValues: Map[DxName, V] =
      WdlUtils.getOutputClosure(outputParamMap).collect {
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
      val paramStrs = outputParams
        .map {
          case (dxName, (wdlType, expr)) =>
            s"${TypeUtils.prettyFormatType(wdlType)} ${dxName} = ${TypeUtils.prettyFormatExpr(expr)}"
        }
        .mkString("\n  ")
      logger.trace(
          s"""|input values: ${inputWdlValues}
              |output closure values: ${outputWdlValues} 
              |outputs parameters:
              |  ${paramStrs}
              |""".stripMargin
      )
    }
    // evaluate
    val env = WdlValueBindings((inputWdlValues ++ outputWdlValues).map {
      case (dxName, v) => dxName.decoded -> v
    })
    val evaluatedOutputValues = evaluator
      .applyAll(outputParams.map {
        case (dxName, te) => (dxName.decoded, te)
      }, env)
      .toMap
      .map {
        case (name, value) => WdlDxName.fromDecodedName(name) -> value
      }
    // convert back to IR
    val irOutputs: Map[DxName, (Type, Value)] = evaluatedOutputValues.map {
      case (name, value) =>
        val wdlType = outputParamMap(name)._1
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

  private def evaluateExpression(expr: TAT.Expr, wdlType: T, env: Map[DxName, (T, V)]): V = {
    evaluator.applyExprAndCoerce(expr, wdlType, Eval.createBindingsFromEnv(env.map {
      case (dxName, tv) => dxName.decoded -> tv
    }))
  }

  private def getBlockOutputs(elements: Vector[TAT.WorkflowElement]): Map[DxName, T] = {
    val (_, outputs) = WdlUtils.getClosureInputsAndOutputs(elements, withField = false)
    outputs.map {
      case (dxName, (wdlType, _)) => dxName -> wdlType
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
      env: Map[DxName, (T, V)]
  ): Map[DxName, (T, V)] = {
    elements.foldLeft(Map.empty[DxName, (T, V)]) {
      case (accu, TAT.PrivateVariable(name, wdlType, expr)) =>
        val value = evaluateExpression(expr, wdlType, accu ++ env)
        accu + (WdlDxName.fromSourceName(name) -> (wdlType, value))
      case (accu, TAT.Conditional(expr, body)) =>
        // evaluate the condition
        val results: Map[DxName, (T, V)] =
          evaluateExpression(expr, expr.wdlType, accu ++ env) match {
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
      case (accu, TAT.Scatter(id, expr, body)) =>
        val collection: Vector[V] =
          evaluateExpression(expr, expr.wdlType, accu ++ env) match {
            case V_Array(array) => array
            case other =>
              throw new AppInternalException(s"Unexpected class ${other.getClass}, ${other}")
          }
        val outputTypes: Map[DxName, T] = getBlockOutputs(body)
        val outputNames: Vector[DxName] = outputTypes.keys.toVector
        val results = if (collection.isEmpty) {
          // empty scatter - return an empty array for each name
          outputNames.map { name =>
            val arrayType = T_Array(outputTypes(name), nonEmpty = false)
            name -> (arrayType, V_Array())
          }
        } else {
          // iterate on the collection, evaluate the body N times, transpose the results into M
          // vectors of N items
          val outputValues: Vector[Vector[V]] =
            collection.map { v =>
              val envInner = accu ++ env + (WdlDxName.fromSourceName(id) -> (expr.wdlType, v))
              val bodyValues = evaluateWorkflowElementVariables(body, envInner)
              outputNames.map(bodyValues(_)._2)
            }.transpose
          // Add the WDL array type to each vector
          outputNames.zip(outputValues).map {
            case (name, values) =>
              val arrayType = T_Array(outputTypes(name), nonEmpty = false)
              val arrayValue = V_Array(values)
              name -> (arrayType, arrayValue)
          }
        }
        accu ++ results
      case (_, other) =>
        throw new Exception(s"type ${other.getClass} while evaluating expressions")
    }
  }

  object WdlBlockContext {
    private[wdl] def evaluateCallInputs(
        call: TAT.Call,
        env: Map[DxName, (T, V)] = Map.empty
    ): Map[DxName, (T, V)] = {
      call.callee.input.flatMap {
        case (name, (wdlType, optional)) if call.inputs.contains(name) =>
          val dxName = WdlDxName.fromSourceName(name)
          val optType = if (optional) {
            TypeUtils.ensureOptional(wdlType)
          } else {
            TypeUtils.unwrapOptional(wdlType)
          }
          val value =
            try {
              evaluateExpression(call.inputs(name), optType, env)
            } catch {
              case ee: EvalException =>
                call.inputs(name) match {
                  case _: TAT.ValueNull if !TypeUtils.isOptional(optType) =>
                    throw new Exception(
                        s"""missing/null value for non-optional input ${name} to call 
                           |${call.fullyQualifiedName}""".stripMargin.replaceAll("\n", " "),
                        ee
                    )
                  case _ =>
                    throw new Exception(
                        s"""Error evaluating input ${name} value ${call.inputs(name)} to call 
                           |${call.fullyQualifiedName}""".stripMargin.replaceAll("\n", " "),
                        ee
                    )
                }
            }
          Some(dxName -> (optType, value))
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

  case class WdlBlockContext(block: WdlBlock,
                             inputEnv: Map[DxName, (T, V)],
                             prereqEnv: Map[DxName, (T, V)])
      extends BlockContext {
    private def call: TAT.Call = block.call
    private lazy val wdlEnv = inputEnv ++ prereqEnv
    private lazy val prereqIrEnv = WdlUtils.toIR(prereqEnv)
    override lazy val env: Map[DxName, (Type, Value)] = WdlUtils.toIR(inputEnv) ++ prereqIrEnv

    private def evaluateCallInputs(
        extraEnv: Map[DxName, (T, V)] = Map.empty
    ): Map[DxName, (T, V)] = {
      WdlBlockContext.evaluateCallInputs(call, wdlEnv ++ extraEnv)
    }

    private def launchCall(
        callInputs: Map[DxName, (T, V)],
        folder: Option[String],
        nameDetail: Option[String] = None
    ): (DxExecution, ExecutableLink, String) = {
      logger.traceLimited(
          s"""|call = ${call}
              |callInputs = ${callInputs}
              |""".stripMargin,
          minLevel = TraceLevel.VVerbose
      )
      val callInputsIR = WdlUtils.toIR(callInputs)
      val callee = tasks.get(call.callee.name)
      val (instanceType, isNative) = callee
        .map { task =>
          val runtimeIsNative = task.runtime.exists(runtime => runtime.kvs.contains("dx_app"))
          val hintsDx = task.hints.flatMap { hints =>
            hints.kvs.get("dnanexus").flatMap {
              case TAT.MetaValueObject(fields) => Some(fields)
              case _                           => None
            }
          }
          val hintsIsNative = hintsDx.exists(_.contains("app"))

          val isNative = runtimeIsNative || hintsIsNative || task.meta.exists { meta =>
            meta.kvs.get("type") match {
              case Some(TAT.MetaValueString("native", _)) => true
              case _                                      => false
            }
          }
          val hasResources = task.runtime.exists(_.kvs.size > (if (runtimeIsNative) 1 else 0)) ||
            hintsDx.exists(_.size > (if (hintsIsNative) 1 else 0))
          val instanceType = if (hasResources) {
            val callIO = TaskInputOutput(task, logger)
            val inputWdlValues: Map[DxName, V] = callInputsIR.collect {
              case (dxName, (t, v))
                  if !dxName.suffix.exists(_.endsWith(Constants.FlatFilesSuffix)) =>
                val wdlType = WdlUtils.fromIRType(t, wdlTypeAliases)
                dxName -> WdlUtils.fromIRValue(v, wdlType, dxName.decoded)
            }
            // add default values for any missing inputs
            val callInputs = callIO.inputsFromValues(inputWdlValues.map {
              case (dxName, v) => dxName.decoded -> v
            }, evaluator, ignoreDefaultEvalError = false)
            val runtime = Runtime(versionSupport.version,
                                  task.runtime,
                                  task.hints,
                                  evaluator,
                                  ctx = Some(callInputs))
            try {
              val request = runtime.parseInstanceType
              if (request.isEmpty && isNative) {
                // TODO: this will never happen - request always gets filled in with default values.
                //  Currently we're relying on the task wrapper not having a runtime or hints section.
                //  We should change `isEmpty` to something like `isDefault`, which would return true
                //  if the request only consists of default values.
                None
              } else {
                val instanceType = jobMeta.instanceTypeDb.apply(request)
                logger.trace(s"Precalculated instance type for ${task.name}: ${instanceType.name}")
                Some(instanceType)
              }
            } catch {
              case e: Throwable =>
                logger.trace(
                    s"""|Failed to precalculate the instance type for task ${task.name}.
                        |${e}
                        |""".stripMargin
                )
                None
            }
          } else {
            // Try to get instance type from the job details (for frag apps may be populated by ApplicationCompiler)
            val instanceType = thisExecDefaultInstance match {
              case Some(default: String) => jobMeta.instanceTypeDb.selectByName(default)
              case None                  => None
            }
            instanceType
          }
          (instanceType, isNative)
        }
        .getOrElse((None, false))
      // TODO: in the case where instanceType is None, we need to signal to the target that it
      //  needs to try to calculate it's own instance type (unless it's a native app(let)), which
      //  probably means adding an additional input parameter.
      val executableLink = getExecutableLink(call.callee.name)
      val (dxExecution, execName) = launchJob(
          executableLink,
          call.actualName,
          callInputsIR,
          extraManifestOutputs = Option.when(prereqEnv.nonEmpty)(prereqIrEnv),
          nameDetail = nameDetail,
          instanceType = instanceType.map(_.name),
          folder = folder,
          prefixOutputs = true
      )
      logger.trace(
          s"""launched call ${call.actualName} to ${call.callee.name}
             |  executable: ${executableLink}
             |  native?: ${isNative}
             |  execution: ${dxExecution.id}
             |  instance type: ${instanceType.map(_.name)}""".stripMargin
      )
      (dxExecution, executableLink, execName)
    }

    override protected def launchBlockCall(blockIndex: Int): Map[DxName, ParameterLink] = {
      val callInputs = evaluateCallInputs()
      val (dxExecution, executableLink, callName) =
        launchCall(callInputs, Some(blockIndex.toString))
      // TODO: right now, we export all outputs of the executable, which may include intermediate
      //  outputs. Should we filter executable outputs using call.wdlType.output?
      jobMeta.createExecutionOutputLinks(dxExecution, executableLink.outputs, Some(callName))
    }

    private val qualifiedNameRegexp = "(.+)\\.(.+)".r

    private def lookup(name: DxName, env: Map[DxName, (T, V)]): Option[V] = {
      def inner(innerName: DxName): Option[V] = {
        if (env.contains(innerName)) {
          Some(env(innerName)._2)
        } else {
          // TODO: use dxName.popDecodedIdentifier instead of a regexp
          innerName.decoded match {
            case qualifiedNameRegexp(lhs, rhs) =>
              inner(WdlDxName.fromDecodedName(lhs)).map {
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
            case _ => None
          }
        }
      }
      inner(name)
    }

    private def prepareSubworkflowInputs(
        executableLink: ExecutableLink,
        extraEnv: Map[DxName, (T, V)] = Map.empty
    ): Map[DxName, (Type, Value)] = {
      val inputEnv = wdlEnv ++ extraEnv
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
        case (dxName, t) =>
          val value = lookup(dxName, inputEnv)
          logger.traceLimited(s"lookupInEnv(${dxName} = ${value})", minLevel = TraceLevel.VVerbose)
          value.map { value =>
            val wdlType = WdlUtils.fromIRType(t, wdlTypeAliases)
            dxName -> (t, WdlUtils.toIRValue(value, wdlType))
          }
      }
    }

    override protected def launchConditional(): Map[DxName, ParameterLink] = {
      val cond = block.target match {
        case Some(TAT.Conditional(expr, _)) =>
          evaluateExpression(expr, T_Boolean, wdlEnv)
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
          launchBlockCall(block.index)
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
          val (dxExecution, _) = launchJob(
              executableLink,
              executableLink.name,
              callInputs,
              extraManifestOutputs = Option.when(prereqEnv.nonEmpty)(prereqIrEnv),
              folder = Some(block.index.toString)
          )
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
      evaluateExpression(expr, expr.wdlType, wdlEnv) match {
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

    // create a short, easy to read, description for a scatter element.
    private[wdl] def getScatterName(item: V, index: Int): String = {
      def inner(innerItem: V): Option[String] = {
        innerItem match {
          case V_File(path)      => Some(truncate(getFileName(path)))
          case V_Directory(path) => Some(truncate(getFileName(path)))
          case V_Optional(x)     => inner(x)
          case _ if EvalUtils.isPrimitive(innerItem) =>
            Some(truncate(EvalUtils.formatPrimitive(innerItem)))
          case V_Pair(l, r) =>
            val ls = inner(l)
            val rs = inner(r)
            (ls, rs) match {
              case (Some(ls1), Some(rs1)) => Some(truncate(s"(${ls1}, ${rs1})"))
              case _                      => None
            }
          case V_Array(array) =>
            val itemStr =
              WorkflowExecutor.getComplexScatterName(array.iterator.map(inner))
            Some(s"[${itemStr}]")
          case V_Map(members) =>
            val memberStr = WorkflowExecutor.getComplexScatterName(
                members.iterator.map {
                  case (k, v) =>
                    (inner(k), inner(v)) match {
                      case (Some(keyStr), Some(valStr)) => Some(s"${keyStr}: ${valStr}")
                      case _                            => None
                    }
                }
            )
            Some(s"{${memberStr}}")
          case V_Object(members) =>
            val memberStr = WorkflowExecutor.getComplexScatterName(
                members.iterator.map {
                  case (k, v) =>
                    inner(v) match {
                      case Some(valStr) => Some(s"${k}: ${valStr}")
                      case _            => None
                    }
                }
            )
            Some(s"{${memberStr}}")
          case V_Struct(name, members) =>
            val memberStr = WorkflowExecutor.getComplexScatterName(
                members.iterator.map {
                  case (k, v) =>
                    inner(v) match {
                      case Some(valStr) => Some(s"${k}: ${valStr}")
                      case _            => None
                    }
                }
            )
            Some(s"${name} ${memberStr}")
          case _ => None
        }
      }
      inner(item) match {
        case Some(i) => i
        case None    => s"index ${index}"
      }
    }

    private def launchScatterCallJobs(identifier: DxName,
                                      itemType: T,
                                      collection: Vector[V]): Vector[DxExecution] = {
      collection.zipWithIndex.map {
        case (item, index) =>
          val callInputs = evaluateCallInputs(Map(identifier -> (itemType, item)))
          val callNameDetail = getScatterName(item, jobMeta.scatterStart + index)
          val (dxExecution, _, _) =
            launchCall(callInputs,
                       Some((jobMeta.scatterStart + index).toString),
                       Some(callNameDetail))
          dxExecution
      }
    }

    private def launchScatterSubblockJobs(identifier: DxName,
                                          itemType: T,
                                          collection: Vector[V]): Vector[DxExecution] = {
      assert(execLinkInfo.size == 1)
      val executableLink = execLinkInfo.values.head
      collection.zipWithIndex.map {
        case (item, index) =>
          val callInputs =
            prepareSubworkflowInputs(executableLink, Map(identifier -> (itemType, item)))
          val callNameDetail = getScatterName(item, jobMeta.scatterStart + index)
          val (dxExecution, _) = launchJob(
              executableLink,
              executableLink.name,
              callInputs,
              nameDetail = Some(callNameDetail),
              folder = Some((jobMeta.scatterStart + index).toString)
          )
          dxExecution
      }
    }

    override protected def launchScatter(): Map[DxName, ParameterLink] = {
      val (identifier, itemType, collection, next) = block.target match {
        case Some(TAT.Scatter(identifier, expr, _)) =>
          val (collection, next) = evaluateScatterCollection(expr)
          val itemType = expr.wdlType match {
            case T_Array(t, _) => t
            case _ =>
              throw new Exception(s"scatter type ${expr.wdlType} is not an array")
          }
          (WdlDxName.fromSourceName(identifier), itemType, collection, next)
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
      if (childJobs.isEmpty) {
        // no jobs were launched
        // APPS-1473. In WDL it should not matter the structure of an array output if the scatter is applied
        // over an empty array. Here we don't want to check if the task wrapped in the scatter produces arrays.
        createEmptyScatterOutputs(checkDepth = false)
      } else {
        // if there are remaining chunks this calls a continue sub-job,
        // otherwise it calls a collect sub-job
        launchNextScatterChunk(childJobs, next)
      }
    }

    override protected lazy val outputTypes: Map[DxName, Type] = block.outputs.map {
      case WdlBlockOutput(dxName, wdlType, _) => dxName -> WdlUtils.toIRType(wdlType)
    }.toMap

    override protected def getScatterOutputs(
        execOutputs: Vector[Option[Map[DxName, JsValue]]],
        execName: Option[String]
    ): (Map[DxName, (Type, Value)], Option[Map[DxName, (Type, Value)]]) = {
      val outputTypes: Map[DxName, (DxName, Type)] = block.kind match {
        case BlockKind.ScatterOneCall =>
          call.callee.output.map {
            case (name, wdlType) =>
              val dxName = WdlDxName.fromSourceName(name)
              val fqn = dxName.pushDecodedNamespace(call.actualName)
              val irType = WdlUtils.toIRType(wdlType)
              fqn -> (dxName, irType)
          }
        case BlockKind.ScatterComplex =>
          assert(execLinkInfo.size == 1)
          execLinkInfo.values.head.outputs.map {
            case (name, irType) => name -> (name, irType)
          }
        case _ =>
          throw new RuntimeException(s"invalid block ${block}")
      }
      val arrayOutputs = outputTypes.map {
        case (fqn, (name, irType)) =>
          val arrayType = TArray(irType)
          val value = createScatterOutputArray(execOutputs, name, irType, execName)
          fqn -> (arrayType, value)
      }
      (arrayOutputs, Option.when(prereqEnv.nonEmpty)(prereqIrEnv))
    }
  }

  private val compoundNameRegexp = ".+[.\\[].+".r

  override def evaluateBlockInputs(
      jobInputs: Map[DxName, (Type, Value)]
  ): WdlBlockContext = {
    val block: WdlBlock =
      Block.getSubBlockAt(WdlBlock.createBlocks(workflow.body), jobMeta.blockPath)
    // convert the inputs to WDL
    val initEnv = jobInputs.map {
      case (dxName, (irType, irValue)) =>
        val wdlType = WdlUtils.fromIRType(irType, wdlTypeAliases)
        val wdlValue = WdlUtils.fromIRValue(irValue, wdlType, dxName.decoded)
        dxName -> (wdlType, wdlValue)
    }
    val inputEnv = block.inputs.foldLeft(initEnv) {
      case (accu, blockInput: WdlBlockInput) if accu.contains(blockInput.name) =>
        accu
      case (accu, _: ComputedBlockInput) =>
        // this is the scatter variable, or some expression that references the scatter
        // variable - we can ignore it since it will be evaluated during launching of
        // the scatter jobs
        accu
      case (accu, RequiredBlockInput(dxName, wdlType))
          if compoundNameRegexp.matches(dxName.decoded) =>
        // the input name is a compound reference - evaluate it as an identifier
        val expr = versionSupport.parseExpression(dxName.decoded, DefaultBindings(accu.map {
          case (dxName, (wdlType, _)) => dxName.decoded -> wdlType
        }), docSource)
        accu + (dxName -> (wdlType, evaluateExpression(expr, wdlType, accu)))
      case (accu, OptionalBlockInput(dxName, wdlType))
          if compoundNameRegexp.matches(dxName.decoded) =>
        try {
          val expr = versionSupport.parseExpression(dxName.decoded, DefaultBindings(accu.map {
            case (dxName, (wdlType, _)) => dxName.decoded -> wdlType
          }), docSource)
          accu + (dxName -> (wdlType, evaluateExpression(expr, wdlType, accu)))
        } catch {
          case ex: EvalException =>
            val errorMsg = exceptionToString(ex, brief = true)
            logger.trace(
                s"""${errorMsg}; setting the input value to null as it is optional""".stripMargin
            )
            accu + (dxName -> (wdlType, V_Null))
          case other: Exception => throw other
        }
      case (_, RequiredBlockInput(name, _)) =>
        throw new Exception(s"missing required input ${name}")
      case (accu, OverridableBlockInputWithStaticDefault(name, wdlType, defaultValue)) =>
        accu + (name -> (wdlType, defaultValue))
      case (accu, OverridableBlockInputWithDynamicDefault(name, wdlType, defaultExpr)) =>
        val wdlValue = evaluateExpression(defaultExpr, wdlType, accu)
        accu + (name -> (wdlType, wdlValue))
      case (accu, OptionalBlockInput(name, wdlType)) =>
        // the input is missing but it could be optional - add a V_ForcedNull
        accu + (name -> (wdlType, V_ForcedNull))
    }
    val inputEnvNew = inputEnv.filter {
      case (_, (_, V_ForcedNull)) => false
      case _                      => true
    }
    val prereqEnv = evaluateWorkflowElementVariables(block.prerequisites, inputEnvNew)
    WdlBlockContext(block, inputEnvNew, prereqEnv)
  }
}
