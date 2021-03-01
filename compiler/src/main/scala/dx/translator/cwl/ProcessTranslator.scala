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
  ExecutableKindWfOutputs,
  Level,
  LinkInput,
  Parameter,
  ParameterAttribute,
  Stage,
  StageInput,
  Value,
  Workflow
}
import dx.core.languages.cwl.{
  CwlBlock,
  CwlBundle,
  CwlDocumentSource,
  CwlUtils,
  DxHints,
  RequirementEvaluator
}
import dx.cwl.{
  CommandInputParameter,
  CommandLineTool,
  CommandOutputParameter,
  CwlSchema,
  CwlType,
  CwlValue,
  Evaluator,
  EvaluatorContext,
  FileValue,
  Process,
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
import dx.util.{FileSourceResolver, Logger}

import java.nio.file.Paths

case class ProcessTranslator(cwlBundle: CwlBundle,
                             typeAliases: Map[String, CwlSchema],
                             locked: Boolean,
                             defaultRuntimeAttrs: Map[String, Value],
                             reorgAttrs: ReorgSettings,
                             perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                             defaultScatterChunkSize: Int,
                             dxApi: DxApi = DxApi.get,
                             fileResolver: FileSourceResolver = FileSourceResolver.get,
                             logger: Logger = Logger.get) {

  private lazy val cwlDefaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = {
    CwlUtils.fromIRValues(defaultRuntimeAttrs, isInput = true)
  }

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

  private case class CwlToolTranslator(tool: CommandLineTool) {
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
            val (actualType, defaultValue) = cwlEvaluator.evaluate(default, input.types, ctx)
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
      Parameter(name, CwlUtils.toIRType(input.types), irDefaultValue, attrs)
    }

    def translateOutput(output: CommandOutputParameter): Parameter = {
      val name = output.id.get.unqualifiedName.get
      val ctx = CwlUtils.createEvaluatorContext(Runtime.empty)
      val irValue = output.outputBinding
        .flatMap(_.outputEval.flatMap { cwlValue =>
          try {
            val (actualType, actualValue) = cwlEvaluator.evaluate(cwlValue, output.types, ctx)
            val (_, value) = CwlUtils.toIRValue(actualValue, actualType)
            Some(value)
          } catch {
            // the expression cannot be statically evaluated
            case _: Throwable => None
          }
        })
      val attrs = translateParameterAttributes(output, hintParameterAttrs)
      Parameter(name, CwlUtils.toIRType(output.types), irValue, attrs)
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
      val docSource = CwlDocumentSource(
          tool.source
            .map(Paths.get(_))
            .getOrElse(
                throw new Exception(s"no document source for tool ${name}")
            )
      )
      val requirementEvaluator = RequirementEvaluator(
          tool.requirements,
          tool.hints,
          Map.empty,
          DxWorkerPaths.default,
          cwlDefaultRuntimeAttrs
      )
      Application(
          name,
          inputs,
          outputs,
          requirementEvaluator.translateInstanceType,
          requirementEvaluator.translateContainer,
          ExecutableKindApplet,
          docSource,
          translateCallableAttributes(tool, hintCallableAttrs),
          requirementEvaluator.translateApplicationRequirements,
          tags = Set("cwl")
      )
    }
  }

  private case class CwlWorkflowTranslator(wf: CwlWorkflow,
                                           availableDependencies: Map[String, Callable])
      extends WorkflowTranslator(wf.name, availableDependencies, reorgAttrs, logger) {
    // Only the toplevel workflow may be unlocked. This happens
    // only if the user specifically compiles it as "unlocked".
    protected lazy val isLocked: Boolean = {
      cwlBundle.primaryCallable match {
        case Some(wf2: CwlWorkflow) => (wf.name != wf2.name) || locked
        case _                      => true
      }
    }

    private lazy val evaluator: Evaluator = Evaluator.default

    private def createWorkflowInput(input: WorkflowInputParameter): (Parameter, Boolean) = {
      val cwlTypes = input.types
      val irType = CwlUtils.toIRType(cwlTypes)
      val (default, isDynamic) = input.default
        .map { d =>
          try {
            val (_, staticValue) = evaluator.evaluate(d, cwlTypes, EvaluatorContext.empty)
            (Some(staticValue), false)
          } catch {
            case _: Throwable => (None, true)
          }
        }
        .getOrElse((None, false))
      (Parameter(input.name, irType, default.map(CwlUtils.toIRValue(_, cwlTypes)._2), Vector.empty),
       isDynamic)
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
      * @param env the environment
      * @return
      */
    private def translateWfFragment(wfName: String,
                                    block: CwlBlock,
                                    blockPath: Vector[Int],
                                    stepPath: Option[String],
                                    env: CallEnv): (Stage, Vector[Callable]) = {}

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
              val (stage, auxCallables) =
                translateWfFragment(wfName, block, blockPath :+ blockNum, stepPath, beforeEnv)
              val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                case (env, param) =>
                  env.add(param.name, (param, LinkInput(stage.dxStage, param.dxName)))
              }
              (stages :+ (stage, auxCallables), afterEnv)
            }
        }
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
        val irType = CwlUtils.toIRType(param.types)
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

    def translateWorkflowLocked(
        name: String,
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        subBlocks: Vector[CwlBlock],
        level: Level.Value
    ): (Workflow, Vector[Callable], Vector[(Parameter, StageInput)]) = ???

    private def translateTopWorkflowLocked(
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        subBlocks: Vector[CwlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      translateWorkflowLocked(wf.name, inputs, outputs, subBlocks, Level.Top)
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
        inputs.map(input => createWorkflowInput(input)._1)
      val commonStageInputs: Vector[StageInput] = inputs.map(_ => EmptyInput)
      val (commonStg, commonApplet) =
        createCommonApplet(wf.name, commonAppletInputs, commonStageInputs, commonAppletInputs)
      val fauxWfInputs: Vector[LinkedVar] = commonStg.outputs.map { param =>
        val stageInput = LinkInput(commonStg.dxStage, param.dxName)
        (param, stageInput)
      }

      val (allStageInfo, env) =
        createWorkflowStages(wf.name, fauxWfInputs, Vector.empty, subBlocks, locked = false)
      val (stages, auxCallables) = allStageInfo.unzip

      // convert the outputs into an applet+stage
      val (outputStage, outputApplet) = createOutputStage(wf.name, outputs, Vector.empty, env)

      val wfInputs = commonAppletInputs.map(param => (param, EmptyInput))
      val wfOutputs =
        outputStage.outputs.map(param => (param, LinkInput(outputStage.dxStage, param.dxName)))
      val wfSource = CwlWorkflowSource()
      val irwf = Workflow(
          wf.name,
          wfInputs,
          wfOutputs,
          commonStg +: stages :+ outputStage,
          wfSource,
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
        val toolTranslator = CwlToolTranslator(tool)
        Vector(toolTranslator.apply)
      case wf: CwlWorkflow =>
        val wfTranslator = CwlWorkflowTranslator(wf, availableDependencies)
        wfTranslator.apply
      case _ =>
        throw new Exception(s"Process type ${process} not supported")
    }
  }
}
