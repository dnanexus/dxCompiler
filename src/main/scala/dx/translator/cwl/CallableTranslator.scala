package dx.translator.cwl

import java.nio.file.Path

import dx.api.DxApi
import dx.core.ir.{Application, Callable, CallableAttribute, ExecutableKindApplet, Parameter, ParameterAttribute, RuntimeRequirement, Value}
import dx.core.languages.cwl.CwlUtils
import dx.translator.{ReorgSettings, RunSpec}
import dx.cwl.{CommandInputParameter, CommandLineTool, CommandOutputParameter, Process, Workflow}


case class CallableTranslator(cwlBundle: CwlBundle,
                              locked: Boolean,
                              defaultRuntimeAttrs: Map[String, Value],
                              reorgAttrs: ReorgSettings,
                              dxApi: DxApi = DxApi.get) {


  def translateCallable(callable: Process, sourceFile: Path): Callable = {
    val filename = sourceFile.getFileName.toString
    callable match {
      case tool: CommandLineTool =>
        val taskTranslator = CwlToolTranslator.apply(tool, sourceFile)
        taskTranslator.apply(tool, sourceFile)
      case wf: Workflow => throw new NotImplementedError("Workflows are not implemented!")
      case other => throw new NotImplementedError(s"Cwl type ${other.getClass} is not supported!")
    }
  }

  def parseCommandOutputParameter(output: CommandOutputParameter): Parameter = {
    val name = output.id.get.name.get
    val dxType = CwlUtils.toIRType(output.types)

    val attributes: Vector[ParameterAttribute] = CwlUtils.createParameterAttributes(output.doc, output.label)

    Parameter(name, dxType, attributes=attributes)
  }

  def parseCommandInputParameter(input: CommandInputParameter): Parameter = {
    val name = input.id.get.name.get
    val dxType = CwlUtils.toIRType(input.types)
    val defaultValue = Option(CwlUtils.toIRValue(input.default.getOrElse(CwlUtils.getDefaultCWLValue(input.types))))
    val attributes: Vector[ParameterAttribute] = CwlUtils.createParameterAttributes(input.doc, input.label)

    Parameter.apply(name.toString, dxType, defaultValue, attributes)
  }

  private case class CwlToolTranslator(tool: Process, sourceFile: Path) {

    def apply(tool: CommandLineTool, sourceFile: Path): Application = {
      val filename = sourceFile.getFileName.toString
      val inputs = tool.inputs.map(parseCommandInputParameter)
      val outputs = tool.outputs.map(parseCommandOutputParameter)
      val instanceType = RunSpec.DefaultInstanceType // FIXME: get from Extras, if not specified use this.
      val container = RunSpec.NoImage // TODO : Requirements - blocked by https://github.com/common-workflow-lab/cwljava/issues/31
      val executableKind = ExecutableKindApplet
      val documentSource = CwlDocumentSource.apply(sourceFile)
      val attributes: Vector[CallableAttribute] = CwlUtils.createCallableAttributes(tool.label, tool.doc) // TODO : Requirements - blocked by https://github.com/common-workflow-lab/cwljava/issues/31
      val requirements = Vector.empty[RuntimeRequirement] // TODO : Requirements - blocked by https://github.com/common-workflow-lab/cwljava/issues/31
      Application(filename, inputs, outputs, instanceType, container, executableKind, documentSource, attributes, requirements)
    }
  }
}
