package dx.translator.cwl

import dx.api.DxApi
import dx.core.io.DxWorkerPaths
import dx.core.ir.{
  Application,
  Callable,
  CallableAttribute,
  EmptyInput,
  ExecutableKindApplet,
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
  Parameter => CwlParameter,
  Workflow => CwlWorkflow
}
import dx.translator.CallableAttributes.{DescriptionAttribute, TitleAttribute}
import dx.translator.ParameterAttributes.{HelpAttribute, LabelAttribute}
import dx.translator.{DxWorkflowAttrs, ReorgSettings, WorkflowTranslator}
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

    private def createWorkflowStages(
        wfName: String,
        wfInputs: Vector[LinkedVar],
        blockPath: Vector[Int],
        subBlocks: Vector[CwlBlock],
        locked: Boolean
    ): (Vector[(Stage, Vector[Callable])], CallEnv) = {}

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
