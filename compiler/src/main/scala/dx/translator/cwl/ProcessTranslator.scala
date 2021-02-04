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
import dx.cwl.{Parameter => CwlParameter, _}
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

  def translateProcess(process: Process): Vector[Callable] = {
    process match {
      case tool: CommandLineTool =>
        val toolTranslator = CwlToolTranslator(tool)
        Vector(toolTranslator.apply)
      case _ =>
        throw new Exception(s"Process type ${process} not supported")
    }
  }
}
