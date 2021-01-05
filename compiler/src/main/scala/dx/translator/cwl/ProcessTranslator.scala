package dx.translator.cwl

import dx.api.DxApi
import dx.core.io.DxWorkerPaths
import dx.core.ir.{
  Application,
  Callable,
  CallableAttribute,
  ExecutableKindApplet,
  Parameter,
  ParameterAttribute,
  Value
}
import dx.core.languages.cwl.{CwlDocumentSource, CwlUtils, DxHints, RequirementEvaluator}
import dx.cwl._
import dx.translator.CallableAttributes.{DescriptionAttribute, TitleAttribute}
import dx.translator.ParameterAttributes.{HelpAttribute, LabelAttribute}
import dx.translator.{DxWorkflowAttrs, ReorgSettings}
import dx.util.{FileSourceResolver, Logger}

import java.nio.file.Paths

case class ProcessTranslator(typeAliases: Map[String, CwlSchema],
                             locked: Boolean,
                             defaultRuntimeAttrs: Map[String, Value],
                             reorgAttrs: ReorgSettings,
                             perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                             defaultScatterChunkSize: Int,
                             dxApi: DxApi = DxApi.get,
                             fileResolver: FileSourceResolver = FileSourceResolver.get,
                             logger: Logger = Logger.get) {

  private lazy val cwlDefaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = {
    CwlUtils.fromIRValues(defaultRuntimeAttrs)
  }

  private def translateParameterAttributes(
      param: CommandParameter,
      hintParameterAttrs: Map[String, Vector[ParameterAttribute]]
  ): Vector[ParameterAttribute] = {
    val metaAttrs = Vector(param.doc.map(HelpAttribute), param.label.map(LabelAttribute)).flatten
    val hintAttrs = hintParameterAttrs.getOrElse(param.name, Vector.empty)
    metaAttrs ++ hintAttrs
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
    private lazy val cwlEvaluator = Evaluator.create(tool.requirements)
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
      val name = input.id.get.name.get
      val cwlType = input.types match {
        case Vector(cwlType) => cwlType
        case other =>
          throw new Exception(s"Mutliple types are not supported ${other}")
      }
      val (actualType, defaultValue) = input.default match {
        case Some(default) =>
          // input with default
          try {
            val ctx = CwlUtils.createEvauatorContext(Runtime.empty)
            val (actualType, defaultValue) = cwlEvaluator.evaluate(default, Vector(cwlType), ctx)
            (actualType, Some(CwlUtils.toIRValue(defaultValue, actualType)))
          } catch {
            case _: Throwable => (cwlType, None)
          }
        case None => (cwlType, None)
      }
      // a required or optional input
      val irType = CwlUtils.toIRType(actualType)
      val attrs = translateParameterAttributes(input, hintParameterAttrs)
      Parameter(name, irType, defaultValue, attrs)
    }

    def translateOutput(output: CommandOutputParameter): Parameter = {
      val name = output.id.get.name.get
      val cwlType = output.types match {
        case Vector(cwlType) => cwlType
        case other =>
          throw new Exception(s"Mutliple types are not supported ${other}")
      }
      val ctx = CwlUtils.createEvauatorContext(Runtime.empty)
      val (actualType, defaultValue) = output.outputBinding
        .flatMap { binding =>
          binding.outputEval.map { cwlValue =>
            val (actualType, actualValue) = cwlEvaluator.evaluate(cwlValue, Vector(cwlType), ctx)
            (actualType, Some(CwlUtils.toIRValue(actualValue, actualType)))
          }
        }
        .getOrElse((cwlType, None))
      val irType = CwlUtils.toIRType(actualType)
      val attrs = translateParameterAttributes(output, hintParameterAttrs)
      Parameter(name, irType, defaultValue, attrs)
    }

    def apply: Application = {
      val name = tool.id.name.get
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

  private case class CwlWorkflowTranslator(wf: Workflow) {
    def apply: Vector[Callable] = {
      logger.trace(s"Translating workflow ${wf.name}")

    }
  }

  def translateProcess(process: Process): Vector[Callable] = {
    process match {
      case tool: CommandLineTool =>
        val toolTranslator = CwlToolTranslator(tool)
        Vector(toolTranslator.apply)
      case wf: Workflow =>
        val workflowTranslator = CwlWorkflowTranslator(wf)
        workflowTranslator.apply
      case _ =>
        throw new Exception(s"Process type ${process} not supported")
    }
  }
}
