package dx.translator.cwl

import java.nio.file.Path

import dx.api.DxApi
import dx.core.ir.{Application, Callable, CallableAttribute, ExecutableKindApplet, Parameter, ParameterAttribute, RuntimeRequirement, Value}
import dx.core.languages.cwl.CwlUtils
import dx.translator.{ReorgSettings, RunSpec}
import dx.cwl.{BooleanValue, CommandInputParameter, CommandLineTool, CommandOutputParameter, CwlSchema, DockerRequirement, LongValue, Process, SchemaDefRequirement, ToolTimeLimitRequirement, WorkReuseRequirement, Workflow}
import dx.translator.RunSpec.{ContainerImage, IgnoreReuseRequirement, TimeoutRequirement}


case class CallableTranslator(cwlBundle: CwlBundle,
                              locked: Boolean,
                              defaultRuntimeAttrs: Map[String, Value],
                              reorgAttrs: ReorgSettings,
                              dxApi: DxApi = DxApi.get) {


  def translateCallable(callable: Process, sourceFile: Path): Callable = {
    callable match {
      case tool: CommandLineTool =>
        val taskTranslator = CwlToolTranslator(tool, sourceFile)
        taskTranslator.apply
      case _: Workflow => throw new NotImplementedError("Workflows are not implemented!")
      case other => throw new NotImplementedError(s"CWL type ${other.getClass} is not supported!")
    }
  }

  def parseCommandOutputParameter(output: CommandOutputParameter, typeAliases: Map[String, CwlSchema] = Map.empty): Parameter = {
    val name = output.id.get.name.get
    val dxType = CwlUtils.toIRType(output.types)
    val attributes: Vector[ParameterAttribute] = CwlUtils.createParameterAttributes(output.doc, output.label)

    Parameter(name, dxType, attributes=attributes)
  }

  def parseCommandInputParameter(input: CommandInputParameter, typeAliases: Map[String, CwlSchema] = Map.empty): Parameter = {
    val name = input.id.get.name.get
    val dxType = CwlUtils.toIRType(input.types)
    val defaultValue = Option(CwlUtils.toIRValue(input.default.getOrElse(CwlUtils.getDefaultCWLValue(input.types))))
    val attributes: Vector[ParameterAttribute] = CwlUtils.createParameterAttributes(input.doc, input.label)

    Parameter.apply(name.toString, dxType, defaultValue, attributes)
  }

  private case class CwlToolTranslator(tool: CommandLineTool, sourceFile: Path) {

    def apply(): Application = {
      val typeAliases = CwlUtils.getTypeAliases(tool.requirements.collect({case a: SchemaDefRequirement => a}))
      val filename = sourceFile.getFileName.toString
      val inputs = tool.inputs.map(parseCommandInputParameter(_, typeAliases))
      val outputs = tool.outputs.map(parseCommandOutputParameter(_, typeAliases))
      val instanceType = RunSpec.DefaultInstanceType
      var container: ContainerImage = RunSpec.NoImage
      val executableKind = ExecutableKindApplet
      val documentSource = CwlDocumentSource.apply(sourceFile)
      val attributes: Vector[CallableAttribute] = CwlUtils.createCallableAttributes(tool.label, tool.doc)
      var requirements: Vector[RuntimeRequirement] = Vector.empty[RuntimeRequirement]
      for (req <- tool.requirements) {
        req match {
          case WorkReuseRequirement(enable) => requirements = requirements :+ IgnoreReuseRequirement(enable.asInstanceOf[BooleanValue].value)
          case ToolTimeLimitRequirement(timeLimit) => requirements = requirements :+ TimeoutRequirement(minutes = Some((timeLimit.asInstanceOf[LongValue].value / 60.0).ceil.toLong)) // todo? test with an expr and number. Currently not testable - problem in underlying java parser
          case _: DockerRequirement => container = RunSpec.NetworkDockerImage // TODO: can there be a dxfile?
          case other => System.err.println(s"${other} was not implemented yet.") // TODO: do we need tro implement everything?
        }
      }
      Application(filename, inputs, outputs, instanceType, container, executableKind, documentSource, attributes, requirements)
    }
  }
}
