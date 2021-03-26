package dx.core.languages.cwl

import dx.cwl.{
  CommandLineTool,
  CwlSchema,
  ExpressionTool,
  Hint,
  HintUtils,
  Identifier,
  Process,
  Requirement,
  Workflow
}
import org.w3id.cwl.cwl1_2.CWLVersion

import scala.collection.immutable.{SeqMap, TreeSeqMap}

case class CwlBundle(version: CWLVersion,
                     primaryProcess: Process,
                     tools: Map[String, CommandLineTool],
                     expressions: Map[String, ExpressionTool],
                     workflows: Map[String, Workflow],
                     requirements: Map[String, Vector[Requirement]],
                     hints: Map[String, Vector[Hint]],
                     processNames: Set[String]) {

  lazy val typeAliases: Map[String, CwlSchema] =
    HintUtils.getSchemaDefs(primaryProcess.requirements)

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
    tools.values.toVector ++ expressions.values.toVector ++ inner(workflows.values).values.toVector
  }
}

object CwlBundle {
  def getProcesses(
      process: Process,
      tools: Map[Identifier, CommandLineTool] = Map.empty,
      expressions: Map[Identifier, ExpressionTool] = Map.empty,
      workflows: Map[Identifier, Workflow] = Map.empty,
      requirements: Map[Identifier, Vector[Requirement]] = Map.empty,
      hints: Map[Identifier, Vector[Hint]] = Map.empty,
      inheritedRequirements: Vector[Requirement] = Vector.empty,
      inheritedHints: Vector[Hint] = Vector.empty
  ): (Map[Identifier, CommandLineTool],
      Map[Identifier, ExpressionTool],
      Map[Identifier, Workflow],
      Map[Identifier, Vector[Requirement]],
      Map[Identifier, Vector[Hint]]) = {
    if (process.id.isEmpty) {
      throw new Exception(s"missing id for process ${process}")
    }
    val newReqs = if (inheritedRequirements.nonEmpty) {
      requirements + (process.id.get -> inheritedRequirements)
    } else {
      requirements
    }
    val newHints = if (inheritedHints.nonEmpty) {
      hints + (process.id.get -> inheritedHints)
    } else {
      hints
    }
    process match {
      case tool: CommandLineTool if tools.contains(tool.id.get) =>
        (tools, expressions, workflows, newReqs, newHints)
      case tool: CommandLineTool =>
        (tools + (tool.id.get -> tool), expressions, workflows, newReqs, newHints)
      case tool: ExpressionTool if expressions.contains(tool.id.get) =>
        (tools, expressions, workflows, newReqs, newHints)
      case tool: ExpressionTool =>
        (tools, expressions + (tool.id.get -> tool), workflows, newReqs, newHints)
      case wf: Workflow if workflows.contains(wf.id.get) =>
        (tools, expressions, workflows, newReqs, newHints)
      case wf: Workflow =>
        val newWorkflows = workflows + (wf.id.get -> wf)
        wf.steps.foldLeft(tools, expressions, newWorkflows, newReqs, newHints) {
          case ((toolAccu, exprAccu, wfAccu, reqAccu, hintAccu), step) =>
            getProcesses(step.run,
                         toolAccu,
                         exprAccu,
                         wfAccu,
                         reqAccu,
                         hintAccu,
                         inheritedRequirements ++ wf.requirements,
                         inheritedHints ++ wf.hints)
        }
      case _ =>
        throw new Exception(s"unsupported process ${process}")
    }
  }

  def create(process: Process): CwlBundle = {
    val version = process.cwlVersion.getOrElse(
        // due to https://github.com/common-workflow-lab/cwljava/issues/43
        // the version does not get set when parsing packed workflows, so
        // for now we assume v1.2 and hope for the best
        CWLVersion.V1_2
        //throw new Exception(s"top-level process does not have a version ${process}")
    )
    process match {
      case tool: CommandLineTool =>
        CwlBundle(version,
                  tool,
                  Map(tool.name -> tool),
                  Map.empty,
                  Map.empty,
                  Map.empty,
                  Map.empty,
                  Set(tool.name))
      case tool: ExpressionTool =>
        CwlBundle(version,
                  tool,
                  Map.empty,
                  Map(tool.name -> tool),
                  Map.empty,
                  Map.empty,
                  Map.empty,
                  Set(tool.name))
      case wf: Workflow =>
        val (tools, expressions, workflows, requirements, hints) = getProcesses(wf)
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
        val requirementsByName = requirements.map {
          case (id, reqs) => id.name.get -> reqs
        }
        val hintsByName = hints.map {
          case (id, hints) => id.name.get -> hints
        }
        CwlBundle(version,
                  wf,
                  toolsByName,
                  expressionsByName,
                  workflowsByName,
                  requirementsByName,
                  hintsByName,
                  allNames3)
    }
  }
}
