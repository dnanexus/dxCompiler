package dx.translator.cwl

import dx.api.DxApi
import dx.core.ir.{Application, Callable, CallableAttribute, ExecutableKindApplet, Parameter, ParameterAttribute, RuntimeRequirement, Type, Value}
import dx.core.languages.cwl.CwlUtils.{toIRType, toIRValue, getDefaultCWLValue}
import dx.translator.{ReorgSettings, RunSpec}
import dx.cwl.{CommandInputParameter, CommandLineTool, CommandOutputParameter, Process, Workflow}
import dx.translator.ParameterAttributes.{HelpAttribute, LabelAttribute}


case class CallableTranslator(cwlBundle: CwlBundle,
                              locked: Boolean,
                              defaultRuntimeAttrs: Map[String, Value],
                              reorgAttrs: ReorgSettings,
                              dxApi: DxApi = DxApi.get) {


  def translateCallable(callable: Process, filename: String = "Unknown"): Callable = {
    callable match {
      case tool: CommandLineTool =>
        val taskTranslator = CwlToolTranslator.apply(tool, filename)
        taskTranslator.apply(tool, filename)
      case wf: Workflow => throw new NotImplementedError("Workflows are not implemented!")
      case other => throw new NotImplementedError(s"Cwl type ${other.getClass} is not supported!")
    }
  }

  // CommandInputParameter ->  doc: Option[String] , but in documentation string | array<string>
  def parseCommandOutputParameter(output: CommandOutputParameter): Parameter = {
    val name = output.id.get.name.get // FIXME: add errors -  name has to exist
    val dxType = toIRType(output.types) // FIXME : returns only one type, what if there are multiple types?

    val helpAttribute: Option[HelpAttribute] = output.doc match { // FIXME: code duplication.
      case Some(s: String) => Option(HelpAttribute(s))
      //      case Some(s: Array[String]) => Option(HelpAttribute(s.mkString)) // FIXME - according to doc
      case None => Option.empty[HelpAttribute]
    }
    val labelAttribute: Option[LabelAttribute] = output.label match {
      case Some(s: String) => Option(LabelAttribute(s))
      case None => Option.empty[LabelAttribute]
    }
    val attributes: Vector[ParameterAttribute] = Vector(helpAttribute, labelAttribute).collect { case Some(i: ParameterAttribute) => i }

    Parameter(name, dxType, attributes=attributes)
  }

  def parseCommandInputParameter(input: CommandInputParameter): Parameter = {
    val name = input.id.get.name.get // FIXME: add errors -  name has to exist
    val dxType = toIRType(input.types)
    val defaultValue = Option(toIRValue(input.default.getOrElse(getDefaultCWLValue(input.types))))
//    val defaultValue = Option(toIRValue(input.default.get))

    val helpAttribute: Option[HelpAttribute] = input.doc match {
      case Some(s: String) => Option(HelpAttribute(s))
      //      case Some(s: Array[String]) => Option(HelpAttribute(s.mkString)) // FIXME - according to doc
      case None => Option.empty[HelpAttribute]
    }


    val labelAttribute: Option[LabelAttribute] = input.label match {
      case Some(s: String) => Option(LabelAttribute(s))
      case None => Option.empty[LabelAttribute]
    }
    val attributes: Vector[ParameterAttribute] = Vector(helpAttribute, labelAttribute).collect { case Some(i: ParameterAttribute) => i }

    Parameter.apply(name.toString, dxType, defaultValue, attributes)
  }

  // TODO : name should be part of process!
  private case class CwlToolTranslator(tool: Process, name: String) {

    def apply(tool: CommandLineTool, name: String): Application = {
      val inputs = tool.inputs.map(parseCommandInputParameter)
      val outputs = tool.outputs.map(parseCommandOutputParameter)
      val instanceType = RunSpec.DefaultInstanceType // FIXME: get from Extras, if not specified use this.
      val container = RunSpec.NoImage // FIXME : from extras probably
      val executableKind = ExecutableKindApplet
      val documentSource = CwlDocumentSource.apply(tool)
      val attributes = Vector.empty[CallableAttribute]
      val requirements = Vector.empty[RuntimeRequirement]
      Application(name, inputs, outputs, instanceType, container, executableKind, documentSource, attributes, requirements)
    }
  }
}
