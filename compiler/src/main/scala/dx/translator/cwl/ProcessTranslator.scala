package dx.translator.cwl

import dx.api.DxApi
import dx.core.Constants
import dx.core.Constants.{ReorgStatus, ReorgStatusCompleted}
import dx.core.io.DxWorkerPaths
import dx.core.ir.RunSpec.{DefaultInstanceType, NoImage}
import dx.core.ir.Type.TString
import dx.core.ir.Value.VString
import dx.core.ir.{
  Application,
  BlockKind,
  Callable,
  CallableAttribute,
  EmptyInput,
  ExecutableKindApplet,
  ExecutableKindWfCustomReorgOutputs,
  ExecutableKindWfFragment,
  ExecutableKindWfOutputs,
  Level,
  LinkInput,
  Parameter,
  ParameterAttribute,
  SourceCode,
  Stage,
  StageInput,
  Type,
  Value,
  Workflow,
  WorkflowInput
}
import dx.core.languages.cwl.{CwlBlock, CwlBundle, CwlUtils, DxHints, RequirementEvaluator}
import dx.cwl.{
  CommandInputParameter,
  CommandLineTool,
  CommandOutputParameter,
  CwlType,
  CwlValue,
  Evaluator,
  EvaluatorContext,
  FileValue,
  Hint,
  Process,
  Requirement,
  Runtime,
  WorkflowInputParameter,
  WorkflowOutputParameter,
  WorkflowStep,
  WorkflowStepInput,
  Parameter => CwlParameter,
  Workflow => CwlWorkflow
}
import dx.translator.CallableAttributes.{DescriptionAttribute, TitleAttribute}
import dx.translator.ParameterAttributes.{HelpAttribute, LabelAttribute}
import dx.translator.{CustomReorgSettings, DxWorkflowAttrs, ReorgSettings, WorkflowTranslator}
import dx.util.{FileSourceResolver, FileUtils, Logger}

import java.nio.file.{Path, Paths}

case class CwlSourceCode(source: Path, override val targets: Vector[String]) extends SourceCode {
  override val language: String = "cwl"
  override def toString: String = FileUtils.readFileContent(source)
}

case class ProcessTranslator(cwlBundle: CwlBundle,
                             locked: Boolean,
                             defaultRuntimeAttrs: Map[String, Value],
                             reorgAttrs: ReorgSettings,
                             perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                             defaultScatterChunkSize: Int,
                             useManifests: Boolean,
                             dxApi: DxApi = DxApi.get,
                             fileResolver: FileSourceResolver = FileSourceResolver.get,
                             logger: Logger = Logger.get) {

  private lazy val cwlDefaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] =
    CwlUtils.fromIRValues(defaultRuntimeAttrs, isInput = true)

  private def translateParameterAttributes(
      param: CwlParameter,
      hintParameterAttrs: Map[String, Vector[ParameterAttribute]]
  ): Vector[ParameterAttribute] = {
    param.getName
      .map { name =>
        val metaAttrs =
          Vector(param.doc.map(HelpAttribute), param.label.map(LabelAttribute)).flatten
        val hintAttrs = hintParameterAttrs.getOrElse(name, Vector.empty)
        metaAttrs ++ hintAttrs
      }
      .getOrElse(Vector.empty)
  }

  private def translateCallableAttributes(
      process: Process,
      hintCallableAttributes: Map[String, Vector[CallableAttribute]]
  ): Vector[CallableAttribute] = {
    val metaAttrs =
      Vector(process.doc.map(DescriptionAttribute), process.label.map(TitleAttribute)).flatten
    val hintAttrs = hintCallableAttributes.getOrElse(process.name, Vector.empty)
    metaAttrs ++ hintAttrs
  }

  private case class CwlToolTranslator(
      tool: CommandLineTool,
      inheritedRequirements: Vector[Requirement] = Vector.empty,
      inheritedHints: Vector[Hint] = Vector.empty
  ) {
    private lazy val cwlEvaluator = Evaluator.create(tool.requirements, tool.hints)
    private lazy val dxHints = tool.hints.collectFirst {
      case dxHints: DxHints => dxHints
    }
    private lazy val hintParameterAttrs: Map[String, Vector[ParameterAttribute]] = {
      dxHints
        .map(_.getParameterAttributes)
        .getOrElse(Map.empty)
    }
    private lazy val hintCallableAttrs: Map[String, Vector[CallableAttribute]] = {
      dxHints.map(_.getCallableAttributes).getOrElse(Map.empty)
    }

    def translateInput(input: CommandInputParameter): Parameter = {
      val name = input.id.get.unqualifiedName.get
      val irDefaultValue = input.default match {
        case Some(default) =>
          try {
            val ctx = CwlUtils.createEvaluatorContext(Runtime.empty)
            val (actualType, defaultValue) = cwlEvaluator.evaluate(default, input.cwlType, ctx)
            defaultValue match {
              case file: FileValue if !CwlUtils.isDxFile(file) =>
                // cannot specify a local file as default - the default will
                // be resolved at runtime
                None
              case _ =>
                val (_, value) = CwlUtils.toIRValue(defaultValue, actualType)
                Some(value)
            }
          } catch {
            case _: Throwable => None
          }
        case None => None
      }
      val attrs = translateParameterAttributes(input, hintParameterAttrs)
      Parameter(name, CwlUtils.toIRType(input.cwlType), irDefaultValue, attrs)
    }

    def translateOutput(output: CommandOutputParameter): Parameter = {
      val name = output.id.get.unqualifiedName.get
      val ctx = CwlUtils.createEvaluatorContext(Runtime.empty)
      val irValue = output.outputBinding
        .flatMap(_.outputEval.flatMap { cwlValue =>
          try {
            val (actualType, actualValue) = cwlEvaluator.evaluate(cwlValue, output.cwlType, ctx)
            val (_, value) = CwlUtils.toIRValue(actualValue, actualType)
            Some(value)
          } catch {
            // the expression cannot be statically evaluated
            case _: Throwable => None
          }
        })
      val attrs = translateParameterAttributes(output, hintParameterAttrs)
      Parameter(name, CwlUtils.toIRType(output.cwlType), irValue, attrs)
    }

    // TODO: If the source CWL has a DockerRequirement with a dockerLoad
    //  URI like dx://dxCompiler_playground:/glnexus_internal, rewrite it
    //  to dx://project-xxx:file-xxx to avoid a runtime lookup
    def apply: Application = {
      val name = tool.name
      val inputs = tool.inputs.collect {
        case i if i.id.exists(_.name.isDefined) => translateInput(i)
      }
      val outputs = tool.outputs.collect {
        case i if i.id.exists(_.name.isDefined) => translateOutput(i)
      }
      val requirementEvaluator = RequirementEvaluator(
          inheritedRequirements ++ tool.requirements,
          inheritedHints ++ tool.hints,
          Map.empty,
          DxWorkerPaths.default,
          cwlDefaultRuntimeAttrs
      )
      val docSource = tool.source.orElse(cwlBundle.primaryProcess.source) match {
        case Some(path) => Paths.get(path)
        case None =>
          throw new Exception(s"no source code for tool ${name}")
      }
      Application(
          name,
          inputs,
          outputs,
          requirementEvaluator.translateInstanceType,
          requirementEvaluator.translateContainer,
          ExecutableKindApplet,
          CwlSourceCode(docSource, Vector(name)),
          translateCallableAttributes(tool, hintCallableAttrs),
          requirementEvaluator.translateApplicationRequirements,
          tags = Set("cwl")
      )
    }
  }

  private case class CwlWorkflowTranslator(
      wf: CwlWorkflow,
      availableDependencies: Map[String, Callable],
      workflowAttrs: Option[DxWorkflowAttrs],
      inheritedRequirements: Vector[Requirement] = Vector.empty,
      inheritedHints: Vector[Hint] = Vector.empty
  ) extends WorkflowTranslator(wf.name, availableDependencies, reorgAttrs, logger) {
    // Only the toplevel workflow may be unlocked. This happens
    // only if the user specifically compiles it as "unlocked".
    protected lazy val isLocked: Boolean = {
      cwlBundle.primaryProcess match {
        case wf2: CwlWorkflow => (wf.name != wf2.name) || locked
        case _                => true
      }
    }

    private lazy val evaluator: Evaluator = Evaluator.default

    override protected lazy val standAloneWorkflow: SourceCode = {
      val docSource = wf.source.orElse(cwlBundle.primaryProcess.source) match {
        case Some(path) => Paths.get(path)
        case None =>
          throw new Exception(s"no source code for tool ${wf.name}")
      }
      CwlSourceCode(docSource, Vector(wf.name))
    }

    private def createWorkflowInput(input: WorkflowInputParameter): Parameter = {
      val cwlType = input.cwlType
      val irType = CwlUtils.toIRType(cwlType)
      val default = input.default
        .flatMap { d =>
          try {
            val (_, staticValue) = evaluator.evaluate(d, cwlType, EvaluatorContext.empty)
            Some(staticValue)
          } catch {
            case _: Throwable => None
          }
        }
      Parameter(input.name, irType, default.map(CwlUtils.toIRValue(_, cwlType)._2), Vector.empty)
    }

    private def callInputToStageInput(callInput: Option[WorkflowStepInput],
                                      calleeParam: Parameter,
                                      env: CallEnv,
                                      locked: Boolean,
                                      callName: String): StageInput = {
      def lookup(key: String): StageInput = {
        env.get(key) match {
          case Some((_, stageInput)) => stageInput
          case None =>
            throw new Exception(
                s"""|input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callName}>
                    |is missing from the environment. We don't have ${key} in the environment.
                    |""".stripMargin.replaceAll("\n", " ")
            )
        }
      }

      callInput match {
        case None if locked =>
          env.log()
          throw new Exception(
              s"""|input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callName}>
                  |is unspecified. This is illegal in a locked workflow.""".stripMargin
                .replaceAll("\n", " ")
          )
        case None =>
          // the argument may be optional, may have a default, or the callee may not
          // use the argument - defer until runtime
          EmptyInput
        case Some(inp) =>
          try {
            lookup(inp.name)
          } catch {
            case _: Throwable if inp.default.isDefined =>
              // the workflow step has a default that will be evaluated at runtime
              EmptyInput
          }
      }
    }

    private def translateCall(call: WorkflowStep, env: CallEnv, locked: Boolean): Stage = {
      val callInputParams = call.inputs.map { param =>
        param.name -> param
      }.toMap
      // Find the callee
      val calleeName = call.name
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
        callInputToStageInput(callInputParams.get(param.name), param, env, locked, call.name)
      }
      Stage(call.name, getStage(), calleeName, inputs, callee.outputVars)
    }

    /**
      * Builds an applet to evaluate a WDL workflow fragment.
      *
      * @param wfName the workflow name
      * @param block the Block to translate into a WfFragment
      * @param blockPath keeps track of which block this fragment represents;
      *                   a top level block is a number. A sub-block of a top-level
      *                   block is a vector of two numbers, etc.
      * @param stepPath the workflow steps that led to this fragment.
      * @param env the environment
      * @return
      */
    private def translateWfFragment(wfName: String,
                                    block: CwlBlock,
                                    blockPath: Vector[Int],
                                    stepPath: Option[String],
                                    env: CallEnv): (Stage, Callable) = {
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

      // Get the applet inputs from the block inputs, which represent the
      // closure required for the block - all the variables defined earlier
      // that are required to evaluate any expression. Some referenced
      // variables may be undefined because they are optional or defined
      // inside the block - we ignore these.
      val (inputParams, stageInputs) = block.inputNames
        .flatMap(env.lookup)
        .map {
          case (fqn, (param, stageInput)) => (param.copy(name = fqn), stageInput)
        }
        .unzip

      // The fragment runner can only handle a single call. If the
      // block already has exactly one call, then we are good. If
      // it contains a scatter/conditional with several calls,
      // then compile the inner block into a sub-workflow. Also
      // Figure out the name of the callable - we need to link with
      // it when we get to the native phase.
      val (stepName, newStepPath) = block.kind match {
        case BlockKind.ExpressionsOnly =>
          (None, None)
        case BlockKind.CallDirect =>
          throw new Exception(s"a direct call should not reach this stage")
        case BlockKind.CallFragment | BlockKind.ConditionalOneCall =>
          // a block with no nested sub-blocks, and a single call, or
          // a conditional with exactly one call in the sub-block
          (Some(block.target.get.name), None)
        case BlockKind.ScatterOneCall =>
          // a scatter with exactly one call in the sub-block
          val stepName = block.target.get.name
          val newScatterPath = stepPath.map(p => s"${p}.${stepName}").getOrElse(stepName)
          (Some(stepName), Some(newScatterPath))
        case _ =>
          throw new Exception(s"unexpected block ${block.prettyFormat}")
      }

      val scatterChunkSize: Option[Int] = newStepPath.map { sctPath =>
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

      val outputParams: Vector[Parameter] = block.outputs.values.map { param =>
        val irType = CwlUtils.toIRType(param.cwlType)
        Parameter(param.name, irType)
      }.toVector

      // create the type map that will be serialized in the applet's details
      val fqnDictTypes: Map[String, Type] = inputParams.map { param: Parameter =>
        param.dxName -> param.dxType
      }.toMap

      val applet = Application(
          s"${wfName}_frag_${getStageId()}",
          inputParams.toVector,
          outputParams,
          DefaultInstanceType,
          NoImage,
          ExecutableKindWfFragment(stepName.toVector, blockPath, fqnDictTypes, scatterChunkSize),
          standAloneWorkflow
      )

      (Stage(stageName, getStage(), applet.name, stageInputs.toVector, outputParams), applet)
    }

    private def createWorkflowStages(
        wfName: String,
        wfInputs: Vector[LinkedVar],
        blockPath: Vector[Int],
        subBlocks: Vector[CwlBlock],
        stepPath: Option[String],
        locked: Boolean
    ): (Vector[(Stage, Vector[Callable])], CallEnv) = {
      logger.trace(s"Assembling workflow backbone $wfName")

      val inputEnv: CallEnv = CallEnv.fromLinkedVars(wfInputs)

      val logger2 = logger.withIncTraceIndent()
      logger2.trace(s"inputs: ${inputEnv.keys}")

      // link together all the stages into a linear workflow
      val (allStageInfo, stageEnv): (Vector[(Stage, Vector[Callable])], CallEnv) =
        subBlocks.zipWithIndex.foldLeft((Vector.empty[(Stage, Vector[Callable])], inputEnv)) {
          case ((stages, beforeEnv), (block: CwlBlock, blockNum: Int)) =>
            if (block.kind == BlockKind.CallDirect) {
              block.target match {
                case Some(step) if CwlUtils.isSimpleCall(step) =>
                  // The block contains exactly one call, with no extra variables.
                  // All the variables are already in the environment, so there is no
                  // need to do any extra work. Compile directly into a workflow stage.
                  logger2.trace(s"Translating step ${step.name} as stage")
                  val stage = translateCall(step, beforeEnv, locked)
                  // Add bindings for the output variables. This allows later calls to refer
                  // to these results.
                  val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                    case (env, param: Parameter) =>
                      val fqn = s"${step.name}/${param.name}"
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
              val (stage, callable) =
                translateWfFragment(wfName, block, blockPath :+ blockNum, stepPath, beforeEnv)
              val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                case (env, param) =>
                  env.add(param.name, (param, LinkInput(stage.dxStage, param.dxName)))
              }
              (stages :+ (stage, Vector(callable)), afterEnv)
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
                                  outputs: Vector[WorkflowOutputParameter],
                                  blockPath: Vector[Int],
                                  env: CallEnv): (Stage, Application) = {
      val paramNames = outputs.flatMap { param =>
        if (param.outputSource.isEmpty) {
          throw new Exception(
              s"there must be at least one 'outputSource' for parameter ${param.name}"
          )
        }
        param.outputSource.map(source => source.replaceAll("/", Parameter.ComplexValueKey))
      }.toSet
      logger.trace(s"paramNames: ${paramNames}")

      val (applicationInputs, stageInputs) = paramNames.map { name =>
        env.get(name) match {
          case Some((param, stageInput)) => (param, stageInput)
          case None =>
            throw new Exception(s"parameter ${name} missing from CallEnv")
        }
      }.unzip

      // build definitions of the output variables - if the expression can be evaluated,
      // set the values as the parameter's default
      val outputParams: Vector[Parameter] = outputs.map { param =>
        val irType = CwlUtils.toIRType(param.cwlType)
        Parameter(param.name, irType, None)
      }

      // Determine kind of application. If a custom reorg app is used and this is a top-level
      // workflow (custom reorg applet doesn't apply to locked workflows), add an output
      // variable for reorg status.
      val (applicationKind, updatedOutputVars) = reorgAttrs match {
        case CustomReorgSettings(_, _, true) if !isLocked =>
          val updatedOutputVars = outputParams :+ Parameter(
              ReorgStatus,
              TString,
              Some(VString(ReorgStatusCompleted))
          )
          (ExecutableKindWfCustomReorgOutputs, updatedOutputVars)
        case _ =>
          (ExecutableKindWfOutputs(blockPath), outputParams)
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
      * and it is always called for nested workflows regarless of whether the top level
      * is locked.
      * @param wfName workflow name
      * @param inputs formal workflow inputs
      * @param outputs workflow outputs
      * @param blockPath the path to the current (sub)workflow, as a vector of block indices
      * @param subBlocks the sub-blocks of the current block
      * @param stepPath the workflow steps that led to this fragment.
      * @param level the workflow level
      * @return
      */
    def translateWorkflowLocked(
        wfName: String,
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        blockPath: Vector[Int],
        subBlocks: Vector[CwlBlock],
        stepPath: Option[String],
        level: Level.Value
    ): (Workflow, Vector[Callable], Vector[(Parameter, StageInput)]) = {
      val wfInputParams = inputs.map(createWorkflowInput)
      val wfInputLinks: Vector[LinkedVar] = wfInputParams.map(p => (p, WorkflowInput(p)))
      val (backboneInputs, commonStageInfo) = if (useManifests) {
        // If we are using manifests, we need an initial applet to merge multiple
        // manifests into a single manifest.
        val commonStageInputs = wfInputParams.map(p => WorkflowInput(p))
        val inputOutputs: Vector[Parameter] = inputs.map { i =>
          Parameter(i.name, CwlUtils.toIRType(i.cwlType))
        }
        val (commonStage, commonApplet) =
          createCommonApplet(wfName, wfInputParams, commonStageInputs, inputOutputs, blockPath)
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
          stepPath,
          locked = true
      )
      val (stages, auxCallables) = (commonStageInfo ++ backboneStageInfo).unzip

      // We need a common output stage for either of two reasons:
      // 1. we need to build an output manifest
      val useOutputStage = useManifests || {
        // 2. an output is used directly as an input
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
        val inputNames = inputs.map(_.name).toSet
        outputs.exists(out => out.outputSource.exists(inputNames.contains))
      }

      val (wfOutputs, finalStages, finalCallables) = if (useOutputStage) {
        val (outputStage, outputApplet) = createOutputStage(wfName, outputs, blockPath, env)
        val wfOutputs = outputStage.outputs.map { param =>
          (param, LinkInput(outputStage.dxStage, param.dxName))
        }
        (wfOutputs, stages :+ outputStage, auxCallables.flatten :+ outputApplet)
      } else {
        val wfOutputs = outputs.map { out =>
          val outputStage = env.get(out.name).map(_._2).getOrElse {
            val Vector(stepName, paramName) = out.outputSource.head.split("/").toVector
            LinkInput(getStage(Some(stepName)), paramName)
          }
          val irType = CwlUtils.toIRType(out.cwlType)
          (Parameter(out.name, irType), outputStage)
        }
        (wfOutputs, stages, auxCallables.flatten)
      }

      (Workflow(wfName,
                wfInputLinks,
                wfOutputs,
                finalStages,
                standAloneWorkflow,
                locked = true,
                level),
       finalCallables,
       wfOutputs)
    }

    private def translateTopWorkflowLocked(
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        subBlocks: Vector[CwlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      translateWorkflowLocked(wf.name, inputs, outputs, Vector.empty, subBlocks, None, Level.Top)
    }

    private def translateTopWorkflowUnlocked(
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        subBlocks: Vector[CwlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      // Create a special applet+stage for the inputs. This is a substitute for
      // workflow inputs. We now call the workflow inputs, "fauxWfInputs" since
      // they are references to the outputs of this first applet.
      val commonAppletInputs: Vector[Parameter] =
        inputs.map(input => createWorkflowInput(input))
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

      val irwf = Workflow(
          wf.name,
          wfInputs,
          wfOutputs,
          commonStg +: stages :+ outputStage,
          standAloneWorkflow,
          locked = false,
          Level.Top
      )
      (irwf, commonApplet +: auxCallables.flatten :+ outputApplet, wfOutputs)
    }

    def translate: (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      logger.trace(s"Translating workflow ${wf.name}")
      // Create a stage per workflow step - steps with conditional and/or scatter
      // require helper applets.
      val subBlocks = CwlBlock.createBlocks(wf)
      if (isLocked) {
        translateTopWorkflowLocked(wf.inputs, wf.outputs, subBlocks)
      } else {
        translateTopWorkflowUnlocked(wf.inputs, wf.outputs, subBlocks)
      }
    }
  }

  def translateProcess(process: Process,
                       availableDependencies: Map[String, Callable]): Vector[Callable] = {
    process match {
      case tool: CommandLineTool =>
        val toolTranslator =
          CwlToolTranslator(tool,
                            cwlBundle.requirements.getOrElse(tool.name, Vector.empty),
                            cwlBundle.hints.getOrElse(tool.name, Vector.empty))
        Vector(toolTranslator.apply)
      case wf: CwlWorkflow =>
        val wfAttrs = perWorkflowAttrs.get(wf.name)
        val wfTranslator = CwlWorkflowTranslator(
            wf,
            availableDependencies,
            wfAttrs,
            cwlBundle.requirements.getOrElse(wf.name, Vector.empty),
            cwlBundle.hints.getOrElse(wf.name, Vector.empty)
        )
        wfTranslator.apply
      case _ =>
        throw new Exception(s"Process type ${process} not supported")
    }
  }
}
