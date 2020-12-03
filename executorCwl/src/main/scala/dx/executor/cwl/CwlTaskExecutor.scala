package dx.executor.cwl

import dx.cwl.{CommandLineTool, CwlType, CwlValue, Parser}
import dx.core.io.StreamFiles.StreamFiles
import dx.core.ir.{Type, Value}
import dx.core.languages.cwl.{CwlUtils, IrToCwlValueBindings, RequirementEvaluator}
import dx.executor.{FileUploader, JobMeta, TaskExecutor}
import dx.util.TraceLevel

object CwlTaskExecutor {
  def create(jobMeta: JobMeta,
             fileUploader: FileUploader,
             streamFiles: StreamFiles): CwlTaskExecutor = {
    if (!Parser.canParse(jobMeta.sourceCode)) {
      throw new Exception(
          s"""source code does not appear to be a CWL document of a supported version
             |${jobMeta.sourceCode}""".stripMargin
      )
    }
    val tool = Parser.parseString(jobMeta.sourceCode) match {
      case tool: CommandLineTool => tool
      case other =>
        throw new Exception(s"expected CWL document to contain a CommandLineTool, not ${other}")
    }
    CwlTaskExecutor(tool, jobMeta, fileUploader, streamFiles)
  }
}

case class CwlTaskExecutor(tool: CommandLineTool,
                           jobMeta: JobMeta,
                           fileUploader: FileUploader,
                           streamFiles: StreamFiles)
    extends TaskExecutor(jobMeta, fileUploader, streamFiles) {

  private val fileResolver = jobMeta.fileResolver
  private val logger = jobMeta.logger

  override def executorName: String = "dxExecutorCwl"

  private lazy val inputTypes: Map[String, Vector[CwlType]] = {
    tool.inputs.collect {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.get.name.get -> param.types
    }.toMap
  }

  private def cwlInputs: Map[String, CwlValue] = {
    // convert IR to CWL values; discard auxiliary fields
    val inputCwlValues: Map[String, CwlValue] = jobMeta.primaryInputs.map {
      case (name, value) =>
        name -> CwlUtils.fromIRValue(value, inputTypes(name), name)
    }
  }

  private def printInputs(inputs: Map[String, CwlValue]): Unit = {
    if (logger.isVerbose) {
      val inputStr = tool.inputs
        .flatMap {
          case param if param.id.forall(_.name.forall(inputs.contains)) =>
            val name = param.id.get.name.get
            Some(s"${name} -> (${param.types}, ${inputs.get(name)})")
          case other =>
            logger.trace(s"no input for parameter ${other}")
            None
        }
        .mkString("\n")
      logger.traceLimited(s"inputs: ${inputStr}")
    }
  }

  private lazy val defaultRuntimeAttrs: Map[String, CwlValue] = {
    CwlUtils.fromIR(jobMeta.defaultRuntimeAttrs)
  }

  private def getRequiredInstanceType(inputs: Map[String, CwlValue] = cwlInputs): String = {
    logger.traceLimited("calcInstanceType", minLevel = TraceLevel.VVerbose)
    printInputs(inputs)
    val env = evaluate(inputs)
    val evaluator =
      RequirementEvaluator(tool.requirements, defaultRuntimeAttrs ++ env, jobMeta.workerPaths)
    val request = evaluator.parseInstanceType
    logger.traceLimited(s"calcInstanceType $request")
    jobMeta.instanceTypeDb.apply(request).name
  }

  override protected lazy val getRequiredInstanceType: String = getRequiredInstanceType()

  override protected def getInputsWithDefaults: Map[String, (Type, Value)] = ???

  override protected def streamFileForInput(parameterName: String): Boolean = ???

  override protected def getSchemas: Map[String, Type.TSchema] = ???

  override protected def writeCommandScript(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = ???

  override protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = ???

  override protected def outputTypes: Map[String, Type] = ???
}
