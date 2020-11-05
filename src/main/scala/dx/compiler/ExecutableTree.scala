package dx.compiler

import dx.api.{DxExecutable, DxWorkflow, Field}
import dx.core.ir.{Application, Callable, ExecutableKind, ExecutableLink, Workflow}
import dx.util.CodecUtils
import spray.json._

case class CompiledExecutable(callable: Callable,
                              dxExec: DxExecutable,
                              links: Vector[ExecutableLink] = Vector.empty,
                              execTree: Option[JsValue] = None)

case class CompilerResults(primary: Option[CompiledExecutable],
                           executables: Map[String, CompiledExecutable]) {
  def executableIds: Vector[String] = {
    primary match {
      case None      => executables.values.map(_.dxExec.id).toVector
      case Some(obj) => Vector(obj.dxExec.id)
    }
  }
}

/**
  * Describe the workflow in a tree representation
  */
case class ExecutableTree(executableDict: Map[String, CompiledExecutable]) {
  def apply(workflow: Workflow): JsValue = {
    val vec = workflow.stages.map { stage =>
      val calleeRecord = executableDict(stage.calleeName)
      val jsv: JsValue = apply(calleeRecord)
      JsObject("stage_name" -> JsString(stage.description), "callee" -> jsv)
    }
    val stages = JsArray(vec)
    JsObject("name" -> JsString(workflow.name), "kind" -> JsString("workflow"), "stages" -> stages)
  }

  def apply(primary: CompiledExecutable): JsValue = {
    primary.callable match {
      case application: Application if primary.links.isEmpty =>
        JsObject("name" -> JsString(application.name),
                 "id" -> JsString(primary.dxExec.id),
                 "kind" -> JsString(ExecutableKind.toString(application.kind)))
      case application: Application =>
        // applet that calls other applets/workflows at runtime.
        // recursively describe all called elements.
        val links: Vector[JsValue] = primary.links.map { eli =>
          val calleeRecord = executableDict(eli.name)
          apply(calleeRecord)
        }
        JsObject(
            "name" -> JsString(application.name),
            "id" -> JsString(primary.dxExec.id),
            "kind" -> JsString(ExecutableKind.toString(application.kind)),
            "executables" -> JsArray(links)
        )
      case workflow: Workflow =>
        val vec = workflow.stages.map { stage =>
          val calleeRecord = executableDict(stage.calleeName)
          val jsv: JsValue = apply(calleeRecord)
          JsObject("stage_name" -> JsString(stage.description), "callee" -> jsv)
        }
        val stages = JsArray(vec)
        JsObject("name" -> JsString(workflow.name),
                 "id" -> JsString(primary.dxExec.id),
                 "kind" -> JsString("workflow"),
                 "stages" -> stages)
    }
  }
}

object ExecutableTree {
  private val CannotFindExecTree = "Unable to find exec tree from"
  private val Indent = 3
  private val LastElement = s"└${"─" * Indent}"
  private val MidElement = s"├${"─" * Indent}"

  private def getDisplayName(stageDesc: Option[String], name: String): String = {
    stageDesc match {
      case Some(name) => name
      case None       => s"${name}"
    }
  }

  private def generateWholePrefix(prefix: String, isLast: Boolean): String = {
    val commonPrefix = prefix.replace("└", " ").replace("─", " ")
    if (isLast) {
      commonPrefix.replace("├", "│") + LastElement
    } else {
      commonPrefix.replace("├", " ") + MidElement
    }
  }

  private def generateTreeBlock(prefix: String,
                                links: Vector[String],
                                title: String,
                                name: String): String = {
    if (links.nonEmpty) {
      prefix + Console.CYAN + title + name + Console.RESET + "\n" + links
        .mkString("\n")
    } else {
      prefix + Console.CYAN + title + name + Console.RESET
    }
  }

  private def processWorkflow(prefix: String, TreeJS: JsObject): String = {
    TreeJS.getFields("name", "stages") match {
      case Seq(JsString(wfName), JsArray(stages)) => {
        val stageLines = stages.zipWithIndex.map {
          case (stage, index) => {
            val isLast = index == stages.length - 1
            val wholePrefix = generateWholePrefix(prefix, isLast)
            val (stageName, callee) = stage.asJsObject.getFields("stage_name", "callee") match {
              case Seq(JsString(stageName), callee: JsObject) => (stageName, callee)
              case x                                          => throw new Exception(s"something is wrong ${x}")
            }
            prettyPrint(callee, Some(stageName), wholePrefix)
          }
        }
        generateTreeBlock(prefix, stageLines, "Workflow: ", Console.YELLOW + wfName)
      }
    }
  }

  private def processApplets(prefix: String,
                             stageDesc: Option[String],
                             TreeJS: JsObject): String = {
    TreeJS.getFields("name", "id", "kind", "executables") match {
      case Seq(JsString(stageName), JsString(_), JsString(kind)) => {
        val name = getDisplayName(stageDesc, stageName)
        generateTreeBlock(prefix, Vector.empty, s"App ${kind}: ", Console.WHITE + name)

      }
      case Seq(JsString(stageName), JsString(_), JsString(kind), JsArray(executables)) => {
        val links = executables.zipWithIndex.map {
          case (link, index) => {
            val isLast = index == (executables.size - 1)
            val wholePrefix = generateWholePrefix(prefix, isLast)
            prettyPrint(link.asJsObject, None, wholePrefix)
          }
        }
        val name = getDisplayName(stageDesc, stageName)
        generateTreeBlock(prefix, links, s"App ${kind}: ", Console.WHITE + name)
      }
      case _ => throw new Exception(s"Missing id, name or kind in ${TreeJS}.")
    }
  }

  /** Recursivly traverse the exec tree and generate an appropriate name + color based on the node type.
    * The prefix is built up as recursive calls happen. This allows for mainaining the locations of branches
    * in the tree. When a prefix made for a current node, it undergoes a transformation to strip out any
    * extra characters from previous calls. This maintains the indenation level and tree branches.
    *
    * Color scheme:
    *   'types' are in CYAN, types being Workflow, App, Task etc
    *   'names' are in WHITE for all app types, and YELLOW for workflow types
    *
    * Step by step:
    *
    *       TREE                                            LEVEL               prettyPrint calls (approx)                                                         NOTES
    * Workflow: four_levels                                 0. Workflow         0a. prettyPrint(IR.Workflow, None)
    * ├───App Inputs: common                                1. App              1a. prettyPrint(IR.App, Some("common"), 3, "├───")
    * ├───App Fragment: if ((username == "a"))              1. App              1b. prettyPrint(IR.AppFrag, Some("if ((username == "a"))", 3, "├───")              The starting prefix for 6 would be ├───└───, that combo gets fixed to be what you actually see by the replace rules
    * │   └───Workflow: four_levels_block_0                 2. Workflow         2a. prettyPrint(IR.Workflow, None, 3, "│   └───")
    * │       ├───App Task: c1                              3. App              3a. prettyPrint(IR.App, None, 3, "│       ├───")                                   The prefix in the step before this would have looked like "│   └───├───"
    * │       └───App Task: c2                              3. App              3b. pretyyPrint(IR.App, None, 3, "│       └───")
    * ├───App Fragment: scatter (i in [1, 4, 9])            1. App              1c. prettyPrint(IR.AppFrag, Some("scatter (i in [1, 4, 9])", 3, "├───")
    * │   └───App Fragment: four_levels_frag_4              2. App              2b. prettyPrint(IR.AppFrag, Some("four_levels_frag_4"), 3, "├───├───")
    * │       └───Workflow: four_levels_block_1_0           3. Workflow          3c. prettyPrint(IR.Workflow, None, 3, "│       └───")
    * │           ├───App Fragment: if ((j == "john"))      4. App              4a. prettyPrint(IR.AppFrag, Some("if ((j == "john"))"), 3, "│           ├───")
    * │           │   └───App Task: concat                  5. App              5a. prettyPrint(IR.App, None, 3, "│           │   └───")                           The prefix that would be 'fixed', into this was "│           ├───└───"
    * │           └───App Fragment: if ((j == "clease"))    4. App              4b. prettyPrint(IR.AppFrag, Some("if ((j == "clease"))"), 3, "│           └───")
    * └───App Outputs: outputs                              1. App              1d. prettyPrint(IR.AppFrag, Some("outputs"), 3, "└───")
    * */
  def prettyPrint(TreeJS: JsObject,
                  stageDesc: Option[String] = None,
                  prefix: String = ""): String = {
    TreeJS.fields.get("kind") match {
      case Some(JsString("workflow")) => processWorkflow(prefix, TreeJS)
      case Some(JsString(_))          => processApplets(prefix, stageDesc, TreeJS)
      case _                          => throw new Exception(s"Missing 'kind' field to be in execTree's entry ${TreeJS}.")
    }
  }

  def fromDxWorkflow(workflow: DxWorkflow): JsValue = {
    val execTree = workflow.describe(Set(Field.Details)).details match {
      case Some(x: JsValue) =>
        x.asJsObject.fields.get("execTree") match {
          case Some(JsString(execString)) => execString
          case _                          => throw new Exception(s"${CannotFindExecTree} for ${workflow.id}")
        }
      case None => throw new Exception(s"${CannotFindExecTree} for ${workflow.id}")
    }
    val TreeJS = CodecUtils.base64DecodeAndGunzip(execTree).parseJson.asJsObject
    JsObject(
        TreeJS.fields + ("id" -> JsString(workflow.id))
    )
  }
}
