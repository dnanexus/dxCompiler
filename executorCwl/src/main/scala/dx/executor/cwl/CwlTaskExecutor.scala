package dx.executor.cwl

import dx.cwl.{CommandLineTool, Parser}
import dx.core.io.StreamFiles.StreamFiles
import dx.core.ir.{Type, Value}
import dx.executor.{FileUploader, JobMeta, TaskExecutor}

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

  override protected def getRequiredInstanceType: String = ???

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
