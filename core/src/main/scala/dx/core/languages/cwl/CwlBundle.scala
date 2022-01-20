package dx.core.languages.cwl

import dx.cwl.{
  CommandLineTool,
  CwlSchema,
  ExpressionTool,
  Hint,
  HintUtils,
  Process,
  Requirement,
  Workflow
}
import org.w3id.cwl.cwl1_2.CWLVersion

import scala.collection.immutable.SeqMap

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
              deps: SeqMap[String, Workflow] = SeqMap.empty): SeqMap[String, Workflow] = {
      wfs.foldLeft(deps) {
        case (accu, wf) =>
          val unsatisfied = wf.steps.map(_.run).collect {
            case wf: Workflow if !accu.contains(wf.simpleName) => wf
          }
          inner(unsatisfied, accu) + (wf.simpleName -> wf)
      }
    }
    // tools have no dependencies so they come first and their ordering doesn't matter
    tools.values.toVector ++ expressions.values.toVector ++ inner(workflows.values).values.toVector
  }
}

object CwlBundle {
  private def getProcesses(
      process: Process,
      tools: Map[String, CommandLineTool] = Map.empty,
      expressions: Map[String, ExpressionTool] = Map.empty,
      workflows: Map[String, Workflow] = Map.empty,
      requirements: Map[String, Vector[Requirement]] = Map.empty,
      hints: Map[String, Vector[Hint]] = Map.empty,
      inheritedRequirements: Vector[Requirement] = Vector.empty,
      inheritedHints: Vector[Hint] = Vector.empty
  ): (Map[String, CommandLineTool],
      Map[String, ExpressionTool],
      Map[String, Workflow],
      Map[String, Vector[Requirement]],
      Map[String, Vector[Hint]]) = {
    if (process.id.isEmpty) {
      throw new Exception(s"missing id for process ${process}")
    }
    val newReqs = if (inheritedRequirements.nonEmpty) {
      requirements + (process.name -> inheritedRequirements)
    } else {
      requirements
    }
    val newHints = if (inheritedHints.nonEmpty) {
      hints + (process.name -> inheritedHints)
    } else {
      hints
    }
    process match {
      case tool: CommandLineTool if tools.contains(tool.name) && tools(tool.name) != tool =>
        throw new Exception(s"two different processes with the same name ${tool.name}")
      case tool: CommandLineTool
          if expressions.contains(tool.name) || workflows.contains(tool.name) =>
        throw new Exception(s"two different processes with the same name ${tool.name}")
      case tool: CommandLineTool if tools.contains(tool.name) =>
        (tools, expressions, workflows, newReqs, newHints)
      case tool: CommandLineTool =>
        (tools + (tool.name -> tool), expressions, workflows, newReqs, newHints)
      case expr: ExpressionTool
          if expressions.contains(expr.name) && expressions(expr.name) != expr =>
        throw new Exception(s"two different processes with the same name ${expr.name}")
      case expr: ExpressionTool if tools.contains(expr.name) || workflows.contains(expr.name) =>
        throw new Exception(s"two different processes with the same name ${expr.name}")
      case expr: ExpressionTool if expressions.contains(expr.name) =>
        (tools, expressions, workflows, newReqs, newHints)
      case expr: ExpressionTool =>
        (tools, expressions + (expr.name -> expr), workflows, newReqs, newHints)
      case wf: Workflow if workflows.contains(wf.name) && workflows(wf.name) != wf =>
        throw new Exception(s"two different processes with the same name ${wf.name}")
      case wf: Workflow if tools.contains(wf.name) || expressions.contains(wf.name) =>
        throw new Exception(s"two different processes with the same name ${wf.name}")
      case wf: Workflow if workflows.contains(wf.name) =>
        (tools, expressions, workflows, newReqs, newHints)
      case wf: Workflow =>
        val newWorkflows = workflows + (wf.name -> wf)
        wf.steps.foldLeft(tools, expressions, newWorkflows, newReqs, newHints) {
          case ((toolAccu, exprAccu, wfAccu, reqAccu, hintAccu), step) =>
            getProcesses(
                step.run.copySimplifyIds(dropNamespace = true,
                                         replacePrefix = (Left(true), None),
                                         simplifyAutoNames = true,
                                         dropCwlExtension = true),
                toolAccu,
                exprAccu,
                wfAccu,
                reqAccu,
                hintAccu,
                inheritedRequirements ++ wf.requirements,
                inheritedHints ++ wf.hints
            )
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
    // get a copy of the process with the ID simplified so that two processes that are identical
    // except for their ID namespace will compare as equal
    process.copySimplifyIds(dropNamespace = true,
                            replacePrefix = (Left(true), None),
                            simplifyAutoNames = true,
                            dropCwlExtension = true) match {
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
        val allNames = tools.keySet ++ expressions.keySet ++ requirements.keySet
        CwlBundle(version, wf, tools, expressions, workflows, requirements, hints, allNames)
    }
  }
}
