package dx.core.languages.cwl

import dx.cwl.{CommandLineTool, ExpressionTool, Identifier, Process, Workflow}
import org.w3id.cwl.cwl1_2.CWLVersion

import scala.collection.immutable.{SeqMap, TreeSeqMap}

case class CwlBundle(version: CWLVersion,
                     primaryCallable: Option[Process],
                     tools: Map[String, CommandLineTool],
                     expressions: Map[String, ExpressionTool],
                     workflows: Map[String, Workflow],
                     processNames: Set[String]) {
  def sortByDependencies: Vector[Process] = {
    def inner(wfs: Iterable[Workflow],
              deps: SeqMap[String, Workflow] = TreeSeqMap.empty): SeqMap[String, Workflow] = {
      wfs.foldLeft(deps) {
        case (accu, wf) =>
          val unsatisfied = wf.steps.map(_.run).collect {
            case wf: Workflow if !accu.contains(wf.name) => wf
          }
          inner(unsatisfied, accu) + (wf.name -> wf)
      }
    }
    // tools have no dependencies so they come first and their ordering doesn't matter
    tools.values.toVector ++ inner(workflows.values).values.toVector
  }
}

object CwlBundle {
  def getProcesses(
      process: Process,
      tools: Map[Identifier, CommandLineTool] = Map.empty,
      expressions: Map[Identifier, ExpressionTool] = Map.empty,
      workflows: Map[Identifier, Workflow] = Map.empty
  )
      : (Map[Identifier, CommandLineTool],
         Map[Identifier, ExpressionTool],
         Map[Identifier, Workflow]) = {
    process match {
      case tool: CommandLineTool if !tools.contains(tool.id) =>
        (tools + (tool.id -> tool), expressions, workflows)
      case tool: ExpressionTool if !expressions.contains(tool.id) =>
        (tools, expressions + (tool.id -> tool), workflows)
      case wf: Workflow if !workflows.contains(wf.id) =>
        wf.steps.foldLeft(tools, expressions, workflows) {
          case ((toolAccu, exprAccu, wfAccu), step) =>
            getProcesses(step.run, toolAccu, exprAccu, wfAccu + (wf.id -> wf))
        }
    }
  }

  def create(process: Process): CwlBundle = {
    val version = process.cwlVersion.getOrElse(
        throw new Exception(s"top-level process does not have a version ${process}")
    )
    process match {
      case tool: CommandLineTool =>
        CwlBundle(version, Some(tool), Map(tool.name -> tool), Map.empty, Map.empty, Set(tool.name))
      case wf: Workflow =>
        val (tools, expressions, workflows) = getProcesses(wf)
        // check that there are no name collisions
        val (toolsByName, allNames) =
          tools.values.foldLeft(Map.empty[String, CommandLineTool], Set.empty[String]) {
            case ((toolAccu, nameAccu), tool) if !nameAccu.contains(tool.name) =>
              (toolAccu + (tool.name -> tool), nameAccu + tool.name)
            case (_, tool) =>
              throw new Exception(s"duplicate name ${tool.name}")
          }
        val (expressionsByName, allNames2) =
          expressions.values.foldLeft(Map.empty[String, ExpressionTool], allNames) {
            case ((exprAccu, nameAccu), expr) if !nameAccu.contains(expr.name) =>
              (exprAccu + (expr.name -> expr), nameAccu + expr.name)
            case (_, expr) =>
              throw new Exception(s"duplicate name ${expr.name}")
          }
        val (workflowsByName, allNames3) =
          workflows.values.foldLeft(Map.empty[String, Workflow], allNames2) {
            case ((wfAccu, nameAccu), wf) if !nameAccu.contains(wf.name) =>
              (wfAccu + (wf.name -> wf), nameAccu + wf.name)
            case (_, wf) =>
              throw new Exception(s"duplicate name ${wf.name}")
          }
        CwlBundle(version, Some(wf), toolsByName, expressionsByName, workflowsByName, allNames3)
    }
  }
}
