package dx.executor.cwl

import dx.core.languages.Language
import dx.core.languages.cwl.{CwlBlock, DxHintSchema}
import dx.cwl.{Parser, Workflow}
import dx.executor.{JobMeta, WorkflowExecutor}
import spray.json.JsString

object CwlWorkflowExecutor {
  def create(jobMeta: JobMeta, separateOutputs: Boolean): CwlWorkflowExecutor = {
    val parser = Parser.create(hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(jobMeta.sourceCode) match {
      case Some((version, "CommandLineTool")) if Language.parse(version) == Language.CwlV1_2 => ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL Workflow document of a supported version
               |${jobMeta.sourceCode}""".stripMargin
        )
    }
    val wfName = jobMeta.getExecutableAttribute("name") match {
      case Some(JsString(name)) => name
      case _                    => throw new Exception("missing executable name")
    }
    val workflow = parser.parseString(jobMeta.sourceCode, name = Some(wfName)) match {
      case tool: Workflow => tool
      case other =>
        throw new Exception(s"expected CWL document to contain a Workflow, not ${other}")
    }
    CwlWorkflowExecutor(workflow, jobMeta, separateOutputs)
  }
}

case class CwlWorkflowExecutor(workflow: Workflow, jobMeta: JobMeta, separateOutputs: Boolean)
    extends WorkflowExecutor[CwlBlock](jobMeta, separateOutputs) {}
