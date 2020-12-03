package dx.executor.cwl

import dx.cwl.{CommandLineTool, CwlOptional, CwlType, CwlValue, NullValue, Parser}
import dx.core.io.StreamFiles.StreamFiles
import dx.core.ir.{Type, Value}
import dx.core.languages.cwl.{CwlEvaluator, CwlUtils, RequirementEvaluator}
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

  private lazy val (inputTypes, defaults): (Map[String, Vector[CwlType]], Map[String, CwlValue]) = {
    val (inputs, defaults) = tool.inputs.collect {
      case param if param.id.forall(_.name.isDefined) =>
        val name = param.id.get.name.get
        (name -> param.types, param.default.map(name -> _))
    }.unzip
    (inputs.toMap, defaults.flatten.toMap)
  }

  private def cwlInputs: Map[String, CwlValue] = {
    val missingTypes = jobMeta.primaryInputs.keySet.diff(inputTypes.keySet)
    if (missingTypes.nonEmpty) {
      throw new Exception(s"no type information given for input(s) ${missingTypes.mkString(",")}")
    }
    // convert IR to CWL values; discard auxiliary fields
    val evaluator = CwlEvaluator(tool.requirements, jobMeta.workerPaths)
    inputTypes.foldLeft(Map.empty[String, CwlValue]) {
      case (env, (name, cwlTypes)) =>
        val cwlValue = jobMeta.primaryInputs.get(name) match {
          case Some(irValue) =>
            CwlUtils.fromIRValue(irValue, cwlTypes, name)
          case None if defaults.contains(name) =>
            evaluator.evaluate(defaults(name), cwlTypes, evaluator.createEvauatorContext(env))
          case None if cwlTypes.exists(CwlOptional.isOptional) =>
            NullValue
          case _ =>
            throw new Exception(s"Missing required input ${name} to tool ${tool.id}")
        }
        env + (name -> cwlValue)
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
