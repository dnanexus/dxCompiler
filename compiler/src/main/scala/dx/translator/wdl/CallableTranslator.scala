package dx.translator.wdl

import dx.api.{DxApi, DxUtils, DxWorkflowStage}
import dx.core.Constants
import dx.core.ir.RunSpec.{DefaultInstanceType, DxFileDockerImage, NoImage}
import dx.translator.{CustomReorgSettings, DefaultReorgSettings, DxWorkflowAttrs, ReorgSettings}
import dx.core.ir._
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.core.Constants.{ReorgStatus, ReorgStatusCompleted}
import dx.core.languages.wdl.{
  CodeGenerator,
  ComputedBlockInput,
  InputKind,
  OptionalBlockInput,
  OverridableBlockInputWithDynamicDefault,
  OverridableBlockInputWithStaticDefault,
  RequiredBlockInput,
  VersionSupport,
  WdlBlock,
  WdlBlockInput,
  WdlBundle,
  WdlDocumentSource,
  WdlUtils,
  WdlWorkflowSource
}
import wdlTools.eval.{DefaultEvalPaths, Eval, EvalException, WdlValueBindings, WdlValues}
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import wdlTools.types.WdlTypes._
import dx.util.{Adjuncts, FileSourceResolver, Logger}
import wdlTools.syntax.Quoting

import scala.annotation.tailrec
import wdlTools.types.TypedAbstractSyntax.ValueString
import wdlTools.types.TypedAbstractSyntax.ValueFile

case class CallableTranslator(wdlBundle: WdlBundle,
                              typeAliases: Map[String, T_Struct],
                              locked: Boolean,
                              defaultRuntimeAttrs: Map[String, Value],
                              reorgAttrs: ReorgSettings,
                              perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                              defaultScatterChunkSize: Int,
                              useManifests: Boolean,
                              instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
                              versionSupport: VersionSupport,
                              dxApi: DxApi = DxApi.get,
                              fileResolver: FileSourceResolver = FileSourceResolver.get,
                              logger: Logger = Logger.get) {

  private lazy val evaluator: Eval =
    Eval(DefaultEvalPaths.empty, Some(wdlBundle.version), Vector.empty, fileResolver, logger)
  private lazy val codegen = CodeGenerator(typeAliases, wdlBundle.version, logger)

  // Return non-local file dependencies of private variables
  private def translateStaticFileDependencies(
      privateVariables: Vector[TAT.PrivateVariable]
  ): Set[String] = {

    // TODO: also consider files nested in arrays, structs

    privateVariables.collect {
      case TAT.PrivateVariable(_, WdlTypes.T_File, ValueFile(value, _)) if value.contains("://") =>
        value
      case TAT.PrivateVariable(_, WdlTypes.T_File, ValueString(value, _, _))
          if value.contains("://") =>
        value
    }.toSet
  }

  private case class WdlTaskTranslator(task: TAT.Task) {
    private lazy val runtime =
      RuntimeTranslator(wdlBundle.version,
                        task.runtime,
                        task.hints,
                        task.meta,
                        defaultRuntimeAttrs,
                        evaluator,
                        dxApi)
    private lazy val adjunctFiles: Vector[Adjuncts.AdjunctFile] =
      wdlBundle.adjunctFiles.getOrElse(task.name, Vector.empty)
    private lazy val meta = ApplicationMetaTranslator(wdlBundle.version, task.meta, adjunctFiles)
    private lazy val parameterMeta =
      ParameterMetaTranslator(wdlBundle.version, task.parameterMeta, task.hints)

    private def translateInput(input: TAT.InputParameter): Parameter = {
      val wdlType = input.wdlType
      val irType = WdlUtils.toIRType(wdlType)
      val attrs = parameterMeta.translateInput(input.name, wdlType)

      input match {
        case TAT.RequiredInputParameter(name, _) => {
          // This is a task "input" parameter of the form:
          //     Int y

          if (isOptional(irType)) {
            throw new Exception(s"Required input ${name} cannot have optional type ${wdlType}")
          }
          Parameter(name, irType, attributes = attrs)
        }
        case TAT.OverridableInputParameterWithDefault(name, _, defaultExpr) =>
          try {
            // This is a task "input" parameter definition of the form:
            //    Int y = 3
            val defaultValue = evaluator.applyConstAndCoerce(defaultExpr, wdlType)
            (wdlType, defaultValue) match {
              case (WdlTypes.T_File | WdlTypes.T_Optional(WdlTypes.T_File), file: WdlValues.V_File)
                  if !WdlUtils.isDxFile(file) =>
                // the default value cannot be specified in the input spec, so we leave it to be
                // evaluated at runtime - e.g. a local file
                Parameter(name, irType, None, attrs)
              case _ =>
                Parameter(name, irType, Some(WdlUtils.toIRValue(defaultValue, wdlType)), attrs)
            }
          } catch {
            // This is a task "input" parameter definition of the form:
            //    Int y = x + 3
            // We treat it as an optional input - the runtime system will
            // evaluate the expression if no value is specified.
            case _: EvalException =>
              val optType = Type.ensureOptional(irType)
              Parameter(name, optType, None, attrs)

          }
        case TAT.OptionalInputParameter(name, _) =>
          val optType = Type.ensureOptional(irType)
          Parameter(name, optType, None, attrs)
      }
    }

    private def translateOutput(output: TAT.OutputParameter, ignoreDefault: Boolean): Parameter = {
      val wdlType = output.wdlType
      val irType = WdlUtils.toIRType(wdlType)
      val defaultValue = if (ignoreDefault || WdlUtils.isPathType(wdlType)) {
        // if the expression is a File, Directory or collection thereof,
        // the paths will be local to the worker so we cannot use them
        // as the default value (since file inputs need to be given as
        // links to dx files)
        None
      } else {
        try {
          val wdlValue = evaluator.applyConstAndCoerce(output.expr, wdlType)
          Some(WdlUtils.toIRValue(wdlValue, wdlType))
        } catch {
          case _: EvalException => None
        }
      }
      val attr = parameterMeta.translateOutput(output.name, wdlType)
      Parameter(output.name, irType, defaultValue, attr)
    }

    def apply: Application = {
      // If the container is stored as a file on the platform, we need to remove
      // the dxURLs in the runtime section, to avoid a runtime lookup. For example:
      //   dx://dxCompiler_playground:/glnexus_internal -> dx://project-xxxx:record-yyyy
      def replaceContainer(runtime: TAT.RuntimeSection,
                           newContainer: String): TAT.RuntimeSection = {
        Set("docker", "container").foreach { key =>
          if (runtime.kvs.contains(key)) {
            return TAT.RuntimeSection(
                runtime.kvs ++ Map(
                    key -> TAT.ValueString(newContainer, T_String, quoting = Quoting.Double)(
                        runtime.kvs(key).loc
                    )
                )
            )(
                runtime.loc
            )
          }
        }
        runtime
      }

      logger.trace(s"Translating task ${task.name}")
      val (kind, isNative) = runtime.translateExecutableKind match {
        case Some(native: ExecutableKindNative) => (native, true)
        case Some(kind)                         => (kind, false)
        case None                               => (ExecutableKindApplet, false)
      }
      val inputs = task.inputs.map(translateInput)
      // the output declarations in a native applet stub have values only because they
      // are required for WDL parsing - they can be safely ignored
      val outputs = task.outputs.map(translateOutput(_, ignoreDefault = isNative))
      val instanceType = runtime.translateInstanceType(instanceTypeSelection)
      val requirements = runtime.translateRequirements
      val staticFileDependencies = translateStaticFileDependencies(task.privateVariables)
      val attributes = meta.translate
      val container = runtime.translateContainer
      val cleanedTask: TAT.Task = container match {
        case DxFileDockerImage(_, dxFile) =>
          val dxURL = DxUtils.dxDataObjectToUri(dxObj = dxFile, includeProject = false)
          task.copy(runtime = task.runtime.map(rt => replaceContainer(rt, dxURL)))(task.loc)
        case _ => task
      }
      val standAloneTask =
        WdlDocumentSource(codegen.createStandAloneTask(cleanedTask), versionSupport)
      Application(
          name = task.name,
          inputs = inputs,
          outputs = outputs,
          instanceType = instanceType,
          container = container,
          kind = kind,
          document = standAloneTask,
          attributes = attributes,
          requirements = requirements,
          staticFileDependencies = staticFileDependencies
      )
    }
  }

  private case class WdlWorkflowTranslator(wf: TAT.Workflow,
                                           availableDependencies: Map[String, Callable],
                                           workflowAttrs: Option[DxWorkflowAttrs]) {
    private lazy val adjunctFiles: Vector[Adjuncts.AdjunctFile] =
      wdlBundle.adjunctFiles.getOrElse(wf.name, Vector.empty)
    private lazy val meta = WorkflowMetaTranslator(wdlBundle.version, wf.meta, adjunctFiles)
    private lazy val parameterMeta = ParameterMetaTranslator(wdlBundle.version, wf.parameterMeta)
    private lazy val standAloneWorkflow = {
      val dependencyNames = WdlUtils.deepFindCalls(wf.body).map(_.unqualifiedName).toSet
      val dependencies =
        availableDependencies.view.filterKeys(dependencyNames.contains).values.toVector
      WdlDocumentSource(codegen.standAloneWorkflow(wf, dependencies), versionSupport)
    }
    // Only the toplevel workflow may be unlocked. This happens
    // only if the user specifically compiles it as "unlocked".
    private lazy val isLocked: Boolean = {
      wdlBundle.primaryCallable match {
        case Some(wf2: TAT.Workflow) =>
          (wf.name != wf2.name) || locked
        case _ =>
          true
      }
    }

    private type LinkedVar = (Parameter, StageInput)

    case class CallEnv(env: Map[String, LinkedVar]) {
      def add(key: String, lvar: LinkedVar): CallEnv = {
        CallEnv(env + (key -> lvar))
      }

      def addAll(items: Vector[(String, LinkedVar)]): CallEnv = {
        CallEnv(env ++ items)
      }

      def contains(key: String): Boolean = env.contains(key)

      def keys: Set[String] = env.keySet

      def get(key: String): Option[LinkedVar] = env.get(key)

      def apply(key: String): LinkedVar = {
        get(key) match {
          case None =>
            log()
            throw new Exception(s"${key} does not exist in the environment.")
          case Some(lvar) => lvar
        }
      }

      def log(): Unit = {
        if (logger.isVerbose) {
          logger.trace("env:")
          val logger2 = logger.withIncTraceIndent()
          stageInputs.map {
            case (name, stageInput) =>
              logger2.trace(s"$name -> ${stageInput}")
          }
        }
      }

      /**
        * Check if the environment has a variable with a binding for
        * a fully-qualified name. For example, if fqn is "A.B.C", then
        * look for "A.B.C", "A.B", or "A", in that order.
        *
        * If the environment has a pair "p", then we want to be able to
        * to return "p" when looking for "p.left" or "p.right".
        *
        * @param fqn fully-qualified name
        * @return
        */
      def lookup(fqn: String): Option[(String, LinkedVar)] = {
        if (env.contains(fqn)) {
          // exact match
          Some(fqn, env(fqn))
        } else {
          // A.B.C --> A.B
          fqn.lastIndexOf(".") match {
            case pos if pos >= 0 => lookup(fqn.substring(0, pos))
            case _               => None
          }
        }
      }

      def stageInputs: Map[String, StageInput] = {
        env.map {
          case (key, (_, stageInput)) => key -> stageInput
        }
      }

      def staticValues: Map[String, (Type, Value)] = {
        env.collect {
          case (key, (Parameter(_, dxType, _, _), StaticInput(value))) =>
            key -> (dxType, value)
          case (key, (_, WorkflowInput(Parameter(_, dxType, Some(value), _)))) =>
            key -> (dxType, value)
        }
      }
    }

    object CallEnv {
      def fromLinkedVars(lvars: Vector[LinkedVar]): CallEnv = {
        CallEnv(lvars.map {
          case (parameter, stageInput) =>
            parameter.name -> (parameter, stageInput)
        }.toMap)
      }
    }

    /**
      * Represents each workflow input with:
      * 1. Parameter, for type declarations
      * 2. StageInput, connecting it to the source of the input
      *
      * It is possible to provide a default value to a workflow input.
      * For example:
      *  workflow w {
      *    Int? x = 3
      *    String s = "glasses"
      *    ...
      *  }
      *
      * @param input the workflow input
      * @return (parameter, isDynamic)
      */
    private def createWorkflowInput(input: WdlBlockInput): (Parameter, Boolean) = {
      val wdlType = input.wdlType
      val irType = WdlUtils.toIRType(wdlType)
      val attr = parameterMeta.translateInput(input.name, input.wdlType)
      input match {
        case RequiredBlockInput(name, _) =>
          if (Type.isOptional(irType)) {
            throw new Exception(s"Required input ${name} cannot have optional type ${wdlType}")
          }
          (Parameter(name, irType, None, attr), false)
        case ComputedBlockInput(name, _) =>
          throw new Exception(s"computed input ${name} cannot be a workflow input")
        case OverridableBlockInputWithStaticDefault(name, _, defaultValue) =>
          (wdlType, defaultValue) match {
            case (WdlTypes.T_File | WdlTypes.T_Optional(WdlTypes.T_File), file: WdlValues.V_File)
                if !WdlUtils.isDxFile(file) =>
              // the default value cannot be specified in the input spec, so we leave it to be
              // evaluated at runtime - e.g. a local file
              (Parameter(name, irType, None, attr), true)
            case _ =>
              (Parameter(name, irType, Some(WdlUtils.toIRValue(defaultValue, wdlType)), attr),
               false)
          }
        case OverridableBlockInputWithDynamicDefault(name, _, _) =>
          // If the default value is an expression that requires evaluation (i.e. not a constant),
          // treat the input as an optional applet input and leave the default value to be calculated
          // at runtime
          val optType = Type.ensureOptional(irType)
          (Parameter(name, optType, None, attr), true)
        case OptionalBlockInput(name, _) =>
          val optType = Type.ensureOptional(irType)
          (Parameter(name, optType, None, attr), false)
      }
    }

    // generate a stage Id, this is a string of the form: 'stage-xxx'
    private val fragNumIter = Iterator.from(0)

    private def getStageId(stageName: Option[String] = None): String = {
      stageName.map(name => s"stage-${name}").getOrElse(s"stage-${fragNumIter.next()}")
    }

    private def getStage(stageName: Option[String] = None): DxWorkflowStage = {
      DxWorkflowStage(getStageId(stageName))
    }

    /**
      * Create a preliminary applet to handle workflow input/outputs. This is
      * used only in the absence of workflow-level inputs/outputs.
      * @param wfName the workflow name
      * @param appletInputs the common applet inputs
      * @param stageInputs the common stage inputs
      * @param outputs the outputs
      * @return
      */
    private def createCommonApplet(wfName: String,
                                   appletInputs: Vector[Parameter],
                                   stageInputs: Vector[StageInput],
                                   outputs: Vector[Parameter],
                                   blockPath: Vector[Int] = Vector.empty): (Stage, Application) = {
      val applet = Application(s"${wfName}_${Constants.CommonStage}",
                               appletInputs,
                               outputs,
                               DefaultInstanceType,
                               NoImage,
                               ExecutableKindWfInputs(blockPath),
                               standAloneWorkflow)
      logger.trace(s"Compiling common applet ${applet.name}")
      val stage = Stage(
          Constants.CommonStage,
          getStage(Some(Constants.CommonStage)),
          applet.name,
          stageInputs,
          outputs
      )
      (stage, applet)
    }

    private def callExprToStageInput(callInputExpr: Option[TAT.Expr],
                                     calleeParam: Parameter,
                                     env: CallEnv,
                                     locked: Boolean,
                                     callFqn: String): StageInput = {

      def lookup(key: String): StageInput = {
        env.get(key) match {
          case Some((_, stageInput)) => stageInput
          case None =>
            throw new Exception(
                s"""|input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callFqn}>
                    |is missing from the environment. We don't have ${key} in the environment.
                    |""".stripMargin.replaceAll("\n", " ")
            )
        }
      }

      callInputExpr match {
        case None if isOptional(calleeParam.dxType) =>
          // optional argument that is not provided
          EmptyInput
        case None if calleeParam.defaultValue.isDefined =>
          // argument that has a default, it can be omitted in the call
          EmptyInput
        case None if locked =>
          env.log()
          throw new Exception(
              s"""|input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callFqn}>
                  |is unspecified. This is illegal in a locked workflow.""".stripMargin
                .replaceAll("\n", " ")
          )
        case None | Some(_: TAT.ValueNone) | Some(_: TAT.ValueNull) =>
          // Perhaps the callee is not going to use the argument, so wait until
          // runtime to determine whether to fail.
          EmptyInput
        case Some(expr) =>
          try {
            // try to evaluate the expression using the constant values from the env
            val paramWdlType = WdlUtils.fromIRType(calleeParam.dxType, typeAliases)
            val bindings = WdlValueBindings(env.staticValues.map {
              case (key, (dxType, value)) =>
                val bindingWdlType = WdlUtils.fromIRType(dxType, typeAliases)
                val bindingValue = WdlUtils.fromIRValue(value, bindingWdlType, key)
                key -> bindingValue
            })
            val value = evaluator.applyExprAndCoerce(expr, paramWdlType, bindings)
            StaticInput(WdlUtils.toIRValue(value, paramWdlType))
          } catch {
            case _: EvalException =>
              // if the expression is an identifier, look it up in the env
              expr match {
                case TAT.ExprIdentifier(id, _) =>
                  lookup(id)
                case TAT.ExprGetName(TAT.ExprIdentifier(id, _), field, _) =>
                  lookup(s"${id}.${field}")
                case _ =>
                  env.log()
                  throw new Exception(
                      s"""|value of input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callFqn}>
                          |cannot be statically evaluated from expression ${expr}.""".stripMargin
                        .replaceAll("\n", " ")
                  )
              }
          }
      }
    }

    /**
      * Compile a call into a stage in an IR Workflow.
      * In a call like:
      *    call lib.native_mk_list as mk_list {
      *      input: a=x, b=5
      *    }
      * it maps callee input `a` to expression `x`. The expressions are
      * trivial because this is a case where we can directly call an applet.
      * @param call the WDL Call object
      * @param env the environment in which the Callee is being called
      * @param locked whether the workflow is locked
      * @return
      */
    private def translateCall(call: TAT.Call, env: CallEnv, locked: Boolean): Stage = {
      // Find the callee
      val calleeName = call.unqualifiedName
      val callee: Callable = availableDependencies.get(calleeName) match {
        case Some(x) => x
        case _ =>
          throw new Exception(
              s"""|Callable ${calleeName} should exist but is missing from the list of known 
                  |tasks/workflows ${availableDependencies.keys}|""".stripMargin
                .replaceAll("\n", " ")
          )
      }
      // Extract the input values/links from the environment
      val inputs: Vector[StageInput] = callee.inputVars.map { param =>
        callExprToStageInput(call.inputs.get(param.name),
                             param,
                             env,
                             locked,
                             call.fullyQualifiedName)
      }
      Stage(call.actualName, getStage(), calleeName, inputs, callee.outputVars)
    }

    /**
      * A block inside a conditional or scatter. If it is simple, we can use a
      * direct call. Otherwise, recursively call into the asssemble-backbone
      * method, and get a locked subworkflow.
      * @param wfName workflow name
      * @param statements statements in the Block
      * @param blockPath block path
      * @param env env
      * @return the generated callable and a Vector of any auxillary callables
      */
    private def translateNestedBlock(wfName: String,
                                     statements: Vector[TAT.WorkflowElement],
                                     blockPath: Vector[Int],
                                     scatterPath: Option[String],
                                     env: CallEnv): (Callable, Vector[Callable]) = {
      val subBlocks = WdlBlock.createBlocks(statements)
      if (subBlocks.size == 1) {
        val block = subBlocks.head
        if (Set(BlockKind.CallDirect, BlockKind.CallWithSubexpressions)
              .contains(block.kind)) {
          throw new RuntimeException("Single call not expected in nested block")
        }
        // At runtime, we will need to execute a workflow fragment. This requires an applet.
        // This is a recursive call, to compile a potentially complex sub-block. It could
        // have many calls generating many applets and subworkflows.
        val (stage, aux) = translateWfFragment(wfName, block, blockPath :+ 0, scatterPath, env)
        val fragName = stage.calleeName
        val main = aux.find(_.name == fragName) match {
          case None    => throw new Exception(s"Could not find $fragName")
          case Some(x) => x
        }
        (main, aux)
      } else {
        // There are several subblocks, we need a subworkflow to string them together.
        // The subworkflow may access variables outside of its scope. For example,
        // stage-0.result, and stage-1.result are inputs to stage-2, that belongs to
        // a subworkflow. Because the subworkflow is locked, we need to make them
        // proper inputs.
        //
        //  |- stage-0.result
        //  |
        //  |- stage-1.result
        //  |
        //  |--- |- stage-2
        //       |
        //       |- stage-3
        //       |
        //       |- stage-4
        //
        val pathStr = blockPath.map(x => x.toString).mkString("_")
        // get the input and output closures of the workflow statements -
        // these are the inputs that come from outside the block, and the
        // outputs that are exposed
        val (statementClosureInputs, statementClosureOutputs) =
          WdlUtils.getClosureInputsAndOutputs(statements, withField = true)
        // create block inputs for the closure inputs
        val allInputs = WdlBlockInput.create(statementClosureInputs)
        val outputs = statementClosureOutputs.values.toVector
        // collect the sub-block inputs that are not workflow inputs or outputs -
        // these are additional inputs from outside the block that need to be
        // supplied as workflow inputs
        val externalNames = (allInputs.map(_.name) ++ outputs.map(_.name)).toSet

        @tailrec
        def containsName(fqn: String): Boolean = {
          if (externalNames.contains(fqn)) {
            // exact match
            true
          } else {
            // A.B.C --> A.B
            fqn.lastIndexOf(".") match {
              case pos if pos >= 0 => containsName(fqn.substring(0, pos))
              case _               => false
            }
          }
        }

        // TODO: will there ever be block inputs that are not included in
        //  statementClosureInputs?
        val closureInputs = subBlocks.flatMap { block =>
          block.inputs.collect {
            case blockInput if !containsName(blockInput.name) =>
              blockInput.name -> (blockInput.wdlType, blockInput.kind)
          }
        }.toMap
        val inputs = allInputs.filter {
          case _: ComputedBlockInput => false
          case _                     => true
        }
        logger.trace(
            s"""|compileNestedBlock
                |    inputs = $inputs
                |    outputs = $outputs
                |    closureInputs = $closureInputs
                |""".stripMargin
        )
        val blockName = s"${wfName}_block_${pathStr}"
        val (subwf, auxCallables, _) = translateWorkflowLocked(
            blockName,
            inputs,
            closureInputs,
            outputs,
            blockPath,
            subBlocks,
            scatterPath,
            Level.Sub
        )
        (subwf, auxCallables)
      }
    }

    /**
      * Builds an applet to evaluate a WDL workflow fragment.
      *
      * @param wfName the workflow name
      * @param block the Block to translate into a WfFragment
      * @param blockPath keeps track of which block this fragment represents;
      *                   a top level block is a number. A sub-block of a top-level
      *                   block is a vector of two numbers, etc.
      * @param env the environment
      * @return
      */
    private def translateWfFragment(wfName: String,
                                    block: WdlBlock,
                                    blockPath: Vector[Int],
                                    scatterPath: Option[String],
                                    env: CallEnv): (Stage, Vector[Callable]) = {
      val stageName = block.getName match {
        case None       => Constants.EvalStage
        case Some(name) => name
      }
      logger.trace(s"Compiling fragment <$stageName> as stage")
      logger
        .withTraceIfContainsKey("GenerateIR")
        .trace(
            s"""|block:
                |${block.prettyFormat}
                |""".stripMargin
        )

      // Figure out the closure required for this block, out of the environment
      // Find the closure of a block, all the variables defined earlier
      // that are required for the calculation.
      //
      // Note: some referenced variables may be undefined. This could be because they are:
      // 1) optional
      // 2) defined -inside- the block
      val closure: Map[String, LinkedVar] = block.inputs.flatMap { i: WdlBlockInput =>
        env.lookup(i.name) match {
          case None               => None
          case Some((name, lVar)) => Some((name, lVar))
        }
      }.toMap

      val inputVars: Vector[Parameter] = closure.map {
        case (fqn, (param: Parameter, _)) =>
          param.copy(name = fqn)
      }.toVector

      // A reversible conversion mul.result --> mul___result. This
      // assumes the '___' symbol is not used anywhere in the original WDL script.
      //
      // This is a simplifying assumption, that is hopefully sufficient. It disallows
      // users from using variables with the ___ character sequence.
      val fqnDictTypes: Map[String, Type] = inputVars.map { param: Parameter =>
        param.dxName -> param.dxType
      }.toMap

      // Figure out the block outputs
      val outputs: Map[String, WdlTypes.T] = block.outputs.map {
        case TAT.OutputParameter(name, wdlType, _) => name -> wdlType
      }.toMap

      // create a Parameter from each block output. The dx:stage
      // will output these Parameter.
      val outputVars = outputs.map {
        case (fqn, wdlType) =>
          val irType = WdlUtils.toIRType(wdlType)
          Parameter(fqn, irType)
      }.toVector

      // The fragment runner can only handle a single call. If the
      // block already has exactly one call, then we are good. If
      // it contains a scatter/conditional with several calls,
      // then compile the inner block into a sub-workflow. Also
      // Figure out the name of the callable - we need to link with
      // it when we get to the native phase.
      logger
        .withTraceIfContainsKey("GenerateIR")
        .trace(s"category : ${block.kind}")

      // complex conditional and scatter blocks may have private variables
      // that need to be added to the environment for use in nested blocks
      lazy val envWithPrivateVars: CallEnv = {
        env.addAll(
            block.prerequisiteVars.map {
              case (name, wdlType) =>
                name -> (Parameter(name, WdlUtils.toIRType(wdlType)), EmptyInput)
            }
        )
      }

      val (innerCall, auxCallables, newScatterPath) =
        block.kind match {
          case BlockKind.ExpressionsOnly =>
            (None, Vector.empty, None)
          case BlockKind.CallDirect =>
            throw new Exception(s"a direct call should not reach this stage")
          case BlockKind.CallWithSubexpressions | BlockKind.CallFragment |
              BlockKind.ConditionalOneCall =>
            // a block with no nested sub-blocks, and a single call, or
            // a conditional with exactly one call in the sub-block
            (Some(block.call.unqualifiedName), Vector.empty, None)
          case BlockKind.ScatterOneCall =>
            // a conditional with exactly one call in the sub-block
            val scatter = block.scatter
            val newScatterPath =
              scatterPath.map(p => s"${p}.${scatter.identifier}").getOrElse(scatter.identifier)
            (Some(block.call.unqualifiedName), Vector.empty, Some(newScatterPath))
          case BlockKind.ConditionalComplex =>
            // a conditional/scatter with multiple calls or other nested elements
            // in the sub-block
            val conditional = block.conditional
            val (callable, aux) =
              translateNestedBlock(wfName,
                                   conditional.body,
                                   blockPath,
                                   scatterPath,
                                   envWithPrivateVars)
            (Some(callable.name), aux :+ callable, None)
          case BlockKind.ScatterComplex =>
            val scatter = block.scatter
            // add the iteration variable to the inner environment
            val varType = scatter.expr.wdlType match {
              case WdlTypes.T_Array(t, _) => WdlUtils.toIRType(t)
              case _ =>
                throw new Exception("scatter doesn't have an array expression")
            }
            val param = Parameter(scatter.identifier, varType)
            val innerEnv = envWithPrivateVars.add(scatter.identifier, (param, EmptyInput))
            val newScatterPath =
              scatterPath.map(p => s"${p}.${scatter.identifier}").getOrElse(scatter.identifier)
            val (callable, aux) =
              translateNestedBlock(wfName, scatter.body, blockPath, Some(newScatterPath), innerEnv)
            (Some(callable.name), aux :+ callable, Some(newScatterPath))
          case _ =>
            throw new Exception(s"unexpected block ${block.prettyFormat}")
        }

      val scatterChunkSize: Option[Int] = newScatterPath.map { sctPath =>
        val scatterChunkSize = workflowAttrs
          .flatMap { wfAttrs =>
            wfAttrs.scatters
              .flatMap {
                _.get(sctPath)
                  .orElse(wfAttrs.scatterDefaults)
                  .flatMap(scatterAttrs => scatterAttrs.chunkSize)
              }
          }
          .getOrElse(defaultScatterChunkSize)
        if (scatterChunkSize > Constants.JobsPerScatterLimit) {
          logger.warning(
              s"The number of jobs per scatter must be between 1-${Constants.JobsPerScatterLimit}"
          )
          Constants.JobsPerScatterLimit
        } else {
          scatterChunkSize
        }
      }

      val applet = Application(
          s"${wfName}_frag_${getStageId()}",
          inputVars,
          outputVars,
          DefaultInstanceType,
          NoImage,
          ExecutableKindWfFragment(innerCall.toVector, blockPath, fqnDictTypes, scatterChunkSize),
          standAloneWorkflow
      )

      val stageInputs: Vector[StageInput] = closure.values.map {
        case (_, stageInput) => stageInput
      }.toVector

      (Stage(stageName, getStage(), applet.name, stageInputs, outputVars), auxCallables :+ applet)
    }

    /**
      * Assembles the backbone of a workflow, having compiled the independent tasks.
      * This is shared between locked and unlocked workflows. Some of the inputs may
      * have default values that are complex expressions, necessitating an initial
      * fragment that performs the evaluation - this only applies to locked workflows,
      * since unlocked workflows always have a "common" applet to handle such expressions.
      *
      * @param wfName workflow name
      * @param wfInputs workflow inputs
      * @param blockPath the path to the current (sub)workflow, as a vector of block indices
      * @param subBlocks the sub-blocks of the current block
      * @param locked whether the workflow is locked
      * @return A tuple (Vector[(Stage, Vector[Callable])], CallEnv), whose first element
      *         is a Vector of stages and the callables included in each stage, and the
      *         second element is a Map of all the input and output variables of the
      *         workflow and its stages.
      */
    private def createWorkflowStages(
        wfName: String,
        wfInputs: Vector[LinkedVar],
        blockPath: Vector[Int],
        subBlocks: Vector[WdlBlock],
        scatterPath: Option[String],
        locked: Boolean
    ): (Vector[(Stage, Vector[Callable])], CallEnv) = {
      logger.trace(s"Assembling workflow backbone $wfName")

      val inputEnv: CallEnv = CallEnv.fromLinkedVars(wfInputs)

      val logger2 = logger.withIncTraceIndent()
      logger2.trace(s"inputs: ${inputEnv.keys}")

      // link together all the stages into a linear workflow
      val (allStageInfo, stageEnv): (Vector[(Stage, Vector[Callable])], CallEnv) =
        subBlocks.zipWithIndex.foldLeft((Vector.empty[(Stage, Vector[Callable])], inputEnv)) {
          case ((stages, beforeEnv), (block: WdlBlock, blockNum: Int)) =>
            if (block.kind == BlockKind.CallDirect) {
              block.target match {
                case Some(call: TAT.Call) =>
                  // The block contains exactly one call, with no extra variables.
                  // All the variables are already in the environment, so there is no
                  // need to do any extra work. Compile directly into a workflow stage.
                  logger2.trace(s"Translating call ${call.actualName} as stage")
                  val stage = translateCall(call, beforeEnv, locked)
                  // Add bindings for the output variables. This allows later calls to refer
                  // to these results.
                  val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                    case (env, param: Parameter) =>
                      val fqn = s"${call.actualName}.${param.name}"
                      val paramFqn = param.copy(name = fqn)
                      env.add(fqn, (paramFqn, LinkInput(stage.dxStage, param.dxName)))
                  }
                  (stages :+ (stage, Vector.empty[Callable]), afterEnv)
                case _ =>
                  throw new Exception(s"invalid DirectCall block ${block}")
              }
            } else {
              // A simple block that requires just one applet, OR
              // a complex block that needs a subworkflow
              val (stage, auxCallables) =
                translateWfFragment(wfName, block, blockPath :+ blockNum, scatterPath, beforeEnv)
              val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                case (env, param) =>
                  env.add(param.name, (param, LinkInput(stage.dxStage, param.dxName)))
              }
              (stages :+ (stage, auxCallables), afterEnv)
            }
        }

      if (logger2.containsKey("GenerateIR")) {
        logger2.trace(s"stages for workflow $wfName = [")
        val logger3 = logger2.withTraceIfContainsKey("GenerateIR", indentInc = 1)
        allStageInfo.foreach {
          case (stage, _) =>
            logger3.trace(
                s"${stage.description}, ${stage.dxStage.id} -> callee=${stage.calleeName}"
            )
        }
        logger2.trace("]")
      }

      (allStageInfo, stageEnv)
    }

    private def createSimpleWorkflowOutput(output: TAT.OutputParameter, env: CallEnv): LinkedVar = {
      val irType = WdlUtils.toIRType(output.wdlType)
      val attr = parameterMeta.translateOutput(output.name, output.wdlType)
      val param = Parameter(output.name, irType, attributes = attr)
      val stageInput: StageInput = if (env.contains(output.name)) {
        env(output.name)._2
      } else {
        try {
          // try to evaluate the output as a constant
          val v = evaluator.applyConstAndCoerce(output.expr, output.wdlType)
          StaticInput(WdlUtils.toIRValue(v, output.wdlType))
        } catch {
          case _: EvalException =>
            output.expr match {
              case TAT.ExprIdentifier(id, _) =>
                // The output is a reference to a previously defined variable
                env(id)._2
              case TAT.ExprGetName(TAT.ExprIdentifier(id2, _), id, _) =>
                // The output is a reference to a previously defined variable
                env(s"$id2.$id")._2
              case _ =>
                // An expression that requires evaluation
                throw new Exception(
                    s"""|Internal error: (${output.expr}) requires evaluation,
                        |which requires constructing an output applet and a stage""".stripMargin
                      .replaceAll("\n", " ")
                )
            }
        }
      }
      (param, stageInput)
    }

    /**
      * Build an applet + workflow stage for evaluating outputs. There are two reasons
      * to be build a special output section:
      * 1. Locked workflow: some of the workflow outputs are expressions.
      *    We need an extra applet+stage to evaluate them.
      * 2. Unlocked workflow: there are no workflow outputs, so we create
      *    them artificially with a separate stage that collects the outputs.
      * @param wfName the workflow name
      * @param outputs the outputs
      * @param env the environment
      * @return
      */
    private def createOutputStage(wfName: String,
                                  outputs: Vector[TAT.OutputParameter],
                                  blockPath: Vector[Int],
                                  env: CallEnv): (Stage, Application) = {
      // split outputs into those that are passed through directly from inputs vs
      // those that require evaluation
      val (outputsToPass, outputsToEval) = outputs.partition(o => env.contains(o.name))
      val outputsToPassEnv = outputsToPass.map(o => o.name -> env(o.name)).toMap
      // create inputs from the closure of the output nodes that need to be evaluate,
      // which includes (recursively) all the variables in the output node expressions
      val outputsToEvalClosureEnv: Map[String, LinkedVar] =
        WdlUtils.getOutputClosure(outputsToEval).keySet.map(name => name -> env(name)).toMap

      val (applicationInputs, stageInputs) =
        (outputsToPassEnv ++ outputsToEvalClosureEnv).values.unzip
      logger.trace(s"inputVars: ${applicationInputs}")

      // build definitions of the output variables - if the expression can be evaluated,
      // set the values as the parameter's default
      val outputVars: Vector[Parameter] = outputs.map {
        case TAT.OutputParameter(name, wdlType, expr) =>
          val value =
            try {
              val v = evaluator.applyConstAndCoerce(expr, wdlType)
              Some(WdlUtils.toIRValue(v, wdlType))
            } catch {
              case _: EvalException => None
            }
          val irType = WdlUtils.toIRType(wdlType)
          val attr = parameterMeta.translateOutput(name, wdlType)
          Parameter(name, irType, value, attr)
      }

      // Determine kind of application. If a custom reorg app is used and this is a top-level
      // workflow (custom reorg applet doesn't apply to locked workflows), add an output
      // variable for reorg status.
      val (applicationKind, updatedOutputVars) = reorgAttrs match {
        case CustomReorgSettings(_, _, true) if !isLocked =>
          val updatedOutputVars = outputVars :+ Parameter(
              ReorgStatus,
              TString,
              Some(VString(ReorgStatusCompleted))
          )
          (ExecutableKindWfCustomReorgOutputs, updatedOutputVars)
        case _ =>
          (ExecutableKindWfOutputs(blockPath), outputVars)
      }
      val application = Application(
          s"${wfName}_${Constants.OutputStage}",
          applicationInputs.toVector,
          updatedOutputVars,
          DefaultInstanceType,
          NoImage,
          applicationKind,
          standAloneWorkflow
      )
      val stage = Stage(
          Constants.OutputStage,
          getStage(Some(Constants.OutputStage)),
          application.name,
          stageInputs.toVector,
          updatedOutputVars
      )
      (stage, application)
    }

    /**
      * Compile a locked workflow. This is called at the top level for locked workflows,
      * and it is always called for nested workflows regardless of whether the top level
      * is locked.
      * @param wfName workflow name
      * @param inputs formal workflow inputs
      * @param closureInputs inputs depended on by `inputs`; mapping of name to (type, optional)
      * @param outputs workflow outputs
      * @param blockPath the path to the current (sub)workflow, as a vector of block indices
      * @param subBlocks the sub-blocks of the current block
      * @param level the workflow level
      * @return
      */
    private def translateWorkflowLocked(
        wfName: String,
        inputs: Vector[WdlBlockInput],
        closureInputs: Map[String, (T, InputKind.InputKind)],
        outputs: Vector[TAT.OutputParameter],
        blockPath: Vector[Int],
        subBlocks: Vector[WdlBlock],
        scatterPath: Option[String],
        level: Level.Level
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      // translate workflow inputs, and also get a Vector of any non-constant
      // expressions that need to be evaluated in the common stage
      val (wfInputParams, dynamicDefaults): (Vector[Parameter], Vector[Boolean]) =
        inputs.map(createWorkflowInput).unzip
      // inputs that are a result of accessing variables in an encompassing
      // WDL workflow.
      val closureInputParams: Vector[Parameter] = closureInputs.map {
        case (name, (wdlType, InputKind.Required)) =>
          // no default value
          val irType = WdlUtils.toIRType(wdlType)
          if (isOptional(irType)) {
            throw new Exception(s"Required input ${name} cannot have optional type ${wdlType}")
          }
          Parameter(name, irType)
        case (name, (wdlType, InputKind.Optional)) =>
          // there is a default value. This input is de facto optional.
          // We change the type of the Parameter and make sure it is optional.
          // no default value
          val irType = WdlUtils.toIRType(wdlType)
          Parameter(name, Type.ensureOptional(irType))
        case (name, (_, InputKind.Computed)) =>
          throw new Exception(s"computed input parameter ${name} not allowed as workflow input")
      }.toVector
      val allWfInputParameters = wfInputParams ++ closureInputParams
      val wfInputLinks: Vector[LinkedVar] = allWfInputParameters.map(p => (p, WorkflowInput(p)))

      val (backboneInputs, commonStageInfo) =
        if (useManifests || dynamicDefaults.exists(identity)) {
          // If we are using manifests, we need an initial applet to merge multiple
          // manifests into a single manifest.
          // If the workflow has inputs that are defined with complex expressions,
          // we need an initial applet to evaluate those.
          val commonAppletInputs = allWfInputParameters
          val commonStageInputs = allWfInputParameters.map(p => WorkflowInput(p))
          val inputOutputs: Vector[Parameter] = inputs.map { i =>
            // TODO: do we need to force the type to be non-optional in the case of
            //  OverridableBlockInputWithDynamicDefault?
            Parameter(i.name, WdlUtils.toIRType(i.wdlType))
          }
          val (commonStage, commonApplet) =
            createCommonApplet(wfName,
                               commonAppletInputs,
                               commonStageInputs,
                               inputOutputs ++ closureInputParams,
                               blockPath)
          val fauxWfInputs: Vector[LinkedVar] = commonStage.outputs.map { param =>
            val link = LinkInput(commonStage.dxStage, param.dxName)
            (param, link)
          }
          (fauxWfInputs, Vector((commonStage, Vector(commonApplet))))
        } else {
          (wfInputLinks, Vector.empty)
        }

      // translate the Block(s) into workflow stages
      val (backboneStageInfo, env) = createWorkflowStages(
          wfName,
          backboneInputs,
          blockPath,
          subBlocks,
          scatterPath,
          locked = true
      )
      val (stages, auxCallables) = (commonStageInfo ++ backboneStageInfo).unzip

      // We need a common output stage for any of three reasons:
      // 1. we need to build an output manifest
      val useOutputStage = useManifests || {
        // 2. any output expressions cannot be resolved without by linking to
        // a workflow input or the output of another stage. Note that a constant
        // value *does* require evaluation because output spec does not allow
        // a default value.
        outputs.exists {
          case TAT.OutputParameter(name, _, _) if env.contains(name) =>
            // The environment has a stage with this output
            false
          case TAT.OutputParameter(_, _, TAT.ExprIdentifier(id, _)) if env.contains(id) =>
            // An identifier that is in scope
            false
          case TAT.OutputParameter(_, _, TAT.ExprGetName(TAT.ExprIdentifier(id2, _), id, _)) =>
            // Access to a field value. For example,
            // c1 is call, and the output section is:
            //  output {
            //     Int? result1 = c1.result
            //     Int? result2 = c2.result
            //  }
            !env.contains(s"$id2.$id")
          case _ => true
        }
      } || {
        // 3. an output is used directly as an input
        // For example, in the small workflow below, 'lane' is used in such a manner.
        //
        // workflow inner {
        //   input {
        //      String lane
        //   }
        //   output {
        //      String blah = lane
        //   }
        // }
        //
        // In locked workflows, it is illegal to access a workflow input directly from
        // a workflow output. It is only allowed to access a stage input/output.
        // In locked workflows, it is illegal to access a workflow input directly from
        // a workflow output. It is only allowed to access a stage input/output.
        val (outputsToPass, outputsToEval) = outputs.partition(o => env.contains(o.name))
        val inputAsOutputNames = outputsToPass.map(_.name).toSet ++
          WdlUtils.getOutputClosure(outputsToEval).keySet
        inputs.map(_.name).toSet.intersect(inputAsOutputNames).nonEmpty
      }

      val (wfOutputs, finalStages, finalCallables) = if (useOutputStage) {
        val (outputStage, outputApplet) = createOutputStage(wfName, outputs, blockPath, env)
        val wfOutputs = outputStage.outputs.map { param =>
          (param, LinkInput(outputStage.dxStage, param.dxName))
        }
        (wfOutputs, stages :+ outputStage, auxCallables.flatten :+ outputApplet)
      } else {
        val wfOutputs = outputs.map(output => createSimpleWorkflowOutput(output, env))
        (wfOutputs, stages, auxCallables.flatten)
      }

      val privateVariables = wf.body.collect {
        case e: TAT.PrivateVariable => e
      }
      val staticFileDependencies = translateStaticFileDependencies(privateVariables)

      (Workflow(
           name = wfName,
           inputs = wfInputLinks,
           outputs = wfOutputs,
           stages = finalStages,
           document = WdlWorkflowSource(wf, versionSupport),
           locked = true,
           level = level,
           attributes = meta.translate,
           staticFileDependencies = staticFileDependencies
       ),
       finalCallables,
       wfOutputs)
    }

    /**
      * Compile a top-level locked workflow.
      */
    private def translateTopWorkflowLocked(
        inputs: Vector[WdlBlockInput],
        outputs: Vector[TAT.OutputParameter],
        subBlocks: Vector[WdlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      translateWorkflowLocked(wf.name,
                              inputs,
                              Map.empty,
                              outputs,
                              Vector.empty,
                              subBlocks,
                              None,
                              Level.Top)
    }

    /**
      * Compile a "regular" (i.e. unlocked) workflow. This function only gets
      * called at the top-level.
      */
    private def translateTopWorkflowUnlocked(
        inputs: Vector[WdlBlockInput],
        outputs: Vector[TAT.OutputParameter],
        subBlocks: Vector[WdlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      // Create a special applet+stage for the inputs. This is a substitute for
      // workflow inputs. We now call the workflow inputs, "fauxWfInputs" since
      // they are references to the outputs of this first applet.
      val commonAppletInputs: Vector[Parameter] =
        inputs.map(input => createWorkflowInput(input)._1)
      val commonStageInputs: Vector[StageInput] = inputs.map(_ => EmptyInput)
      val (commonStg, commonApplet) =
        createCommonApplet(wf.name, commonAppletInputs, commonStageInputs, commonAppletInputs)
      val fauxWfInputs: Vector[LinkedVar] = commonStg.outputs.map { param =>
        val stageInput = LinkInput(commonStg.dxStage, param.dxName)
        (param, stageInput)
      }

      val (allStageInfo, env) =
        createWorkflowStages(wf.name, fauxWfInputs, Vector.empty, subBlocks, None, locked = false)
      val (stages, auxCallables) = allStageInfo.unzip

      // convert the outputs into an applet+stage
      val (outputStage, outputApplet) = createOutputStage(wf.name, outputs, Vector.empty, env)

      val wfInputs = commonAppletInputs.map(param => (param, EmptyInput))
      val wfOutputs =
        outputStage.outputs.map(param => (param, LinkInput(outputStage.dxStage, param.dxName)))
      val wfAttr = meta.translate
      val wfSource = WdlWorkflowSource(wf, versionSupport)

      val privateVariables = wf.body.collect {
        case e: TAT.PrivateVariable => e
      }
      val staticFileDependencies = translateStaticFileDependencies(privateVariables)

      val irwf = Workflow(
          name = wf.name,
          inputs = wfInputs,
          outputs = wfOutputs,
          stages = commonStg +: stages :+ outputStage,
          document = wfSource,
          locked = false,
          level = Level.Top,
          attributes = wfAttr,
          staticFileDependencies = staticFileDependencies
      )
      (irwf, commonApplet +: auxCallables.flatten :+ outputApplet, wfOutputs)
    }

    /**
      * Creates an applet to reorganize the output files. We want to
      * move the intermediate results to a subdirectory.  The applet
      * needs to process all the workflow outputs, to find the files
      * that belong to the final results.
      * @param wfName workflow name
      * @param wfOutputs workflow outputs
      * @return
      */
    private def createReorgStage(wfName: String,
                                 wfOutputs: Vector[LinkedVar]): (Stage, Application) = {
      val applet = Application(
          s"${wfName}_${Constants.ReorgStage}",
          wfOutputs.map(_._1),
          Vector.empty,
          DefaultInstanceType,
          NoImage,
          ExecutableKindWorkflowOutputReorg,
          standAloneWorkflow
      )
      logger.trace(s"Creating output reorganization applet ${applet.name}")
      // Link to the X.y original variables
      val inputs: Vector[StageInput] = wfOutputs.map(_._2)
      val stage =
        Stage(Constants.ReorgStage,
              getStage(Some(Constants.ReorgStage)),
              applet.name,
              inputs,
              Vector.empty[Parameter])
      (stage, applet)
    }

    private def createCustomReorgStage(wfOutputs: Vector[LinkedVar],
                                       appletId: String,
                                       reorgConfigFile: Option[String]): (Stage, Application) = {
      logger.trace(s"Creating custom output reorganization applet ${appletId}")
      val (statusParam, statusStageInput): LinkedVar = wfOutputs.filter {
        case (x, _) => x.name == ReorgStatus
      } match {
        case Vector(lvar) => lvar
        case other =>
          throw new Exception(
              s"Expected exactly one output with name ${ReorgStatus}, found ${other}"
          )
      }
      val configFile: Option[VFile] = reorgConfigFile.map(VFile)
      val appInputs = Vector(
          statusParam,
          Parameter(Constants.ReorgConfig, TFile, configFile)
      )
      val appletKind = ExecutableKindWorkflowCustomReorg(appletId)
      val applet = Application(
          appletId,
          appInputs,
          Vector.empty,
          DefaultInstanceType,
          NoImage,
          appletKind,
          standAloneWorkflow
      )
      // Link to the X.y original variables
      val inputs: Vector[StageInput] = configFile match {
        case Some(x) => Vector(statusStageInput, StaticInput(x))
        case _       => Vector(statusStageInput)
      }
      val stage =
        Stage(Constants.ReorgStage,
              getStage(Some(Constants.ReorgStage)),
              applet.name,
              inputs,
              Vector.empty[Parameter])
      (stage, applet)
    }

    def apply: Vector[Callable] = {
      logger.trace(s"Translating workflow ${wf.name}")
      // Create a stage per workflow body element (variable block, call,
      // scatter block, conditional block)
      val subBlocks = WdlBlock.createBlocks(wf.body)
      // translate workflow inputs/outputs to equivalent classes defined in Block
      val inputs = wf.inputs.map(WdlBlockInput.translate)
      val (irWf, irCallables, irOutputs) =
        if (isLocked) {
          translateTopWorkflowLocked(inputs, wf.outputs, subBlocks)
        } else {
          translateTopWorkflowUnlocked(inputs, wf.outputs, subBlocks)
        }
      // add a reorg applet if necessary
      val (updatedWf, updatedCallables) = reorgAttrs match {
        case DefaultReorgSettings(true) =>
          val (reorgStage, reorgApl) = createReorgStage(wf.name, irOutputs)
          (irWf.copy(stages = irWf.stages :+ reorgStage), irCallables :+ reorgApl)
        case CustomReorgSettings(appUri, reorgConfigFile, true) if !isLocked =>
          val (reorgStage, reorgApl) = createCustomReorgStage(irOutputs, appUri, reorgConfigFile)
          (irWf.copy(stages = irWf.stages :+ reorgStage), irCallables :+ reorgApl)
        case _ =>
          (irWf, irCallables)
      }
      // validate workflow stages
      val allCallableNames = updatedCallables.map(_.name).toSet ++ availableDependencies.keySet
      val invalidStages =
        updatedWf.stages.filterNot(stage => allCallableNames.contains(stage.calleeName))
      if (invalidStages.nonEmpty) {
        val invalidStageDesc = invalidStages
          .map(stage => s"$stage.description, $stage.id.getId -> ${stage.calleeName}")
          .mkString("    ")
        throw new Exception(s"""|One or more stages reference missing tasks:
                                |stages: $invalidStageDesc
                                |callables: $allCallableNames
                                |""".stripMargin)
      }
      updatedCallables :+ updatedWf
    }
  }

  /**
    * Translates a WDL Callable to IR.
    * @param callable the TAT.Callable to translate
    * @param availableDependencies the available Callables upon which `callable` may depend
    * @return ir.Callable
    */
  def translateCallable(callable: TAT.Callable,
                        availableDependencies: Map[String, Callable]): Vector[Callable] = {
    callable match {
      case task: TAT.Task =>
        val taskTranslator = WdlTaskTranslator(task)
        Vector(taskTranslator.apply)
      case wf: TAT.Workflow =>
        val wfAttrs = perWorkflowAttrs.get(wf.name)
        val wfTranslator = WdlWorkflowTranslator(wf, availableDependencies, wfAttrs)
        wfTranslator.apply
    }
  }
}
