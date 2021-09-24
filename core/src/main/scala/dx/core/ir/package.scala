package dx.core.ir

import dx.api.DxWorkflowStage
import dx.core.ir.RunSpec.{ContainerImage, InstanceType}
import dx.util.Enum

trait ParameterAttribute

/**
  * Compile-time representation of an input or output parameter.
  * @param name the fully-qualified parameter name
  * @param dxType parameter data type
  * @param defaultValue default value
  * @param attributes metadata used to encode DNAnexus applet input/output specification
  *                   fields, such as {help, suggestions, patterns}.
  */
case class Parameter(
    name: DxName,
    dxType: Type,
    defaultValue: Option[Value] = None,
    attributes: Vector[ParameterAttribute] = Vector.empty
)

trait CallableAttribute

/**
  * A unified type representing a workflow, app, or applet.
  * This is useful when compiling WDL workflows, because they can
  * call other WDL workflows and applets. This is done using the
  * same syntax.
  */
trait Callable {
  def name: String
  def inputVars: Vector[Parameter]
  def outputVars: Vector[Parameter]
  def attributes: Vector[CallableAttribute]
  def tags: Set[String]
  def properties: Map[String, String]
  def staticFileDependencies: Set[String]
}

/**
  * @param dependencies: the order in which to compile the workflows and tasks.
  *                      The first element in the vector depends on nothing else.
  *                      Each other element (may) depend on any of the previous
  *                      elements.
  */
case class Bundle(primaryCallable: Option[Callable],
                  allCallables: Map[String, Callable],
                  dependencies: Vector[String],
                  typeAliases: Map[String, Type]) {

  override def toString: String = {
    primaryCallable match {
      case Some(c) => s"Bundle[${c.name}]"
      case None if allCallables.size == 1 =>
        s"Bundle[${allCallables.values.head.name}]"
      case None =>
        s"Bundle[${allCallables.values.head.name} and ${allCallables.size - 1} others]"
    }
  }
}

object ExecutableType extends Enum {
  type ExecutableType = Value
  val App, Applet, Workflow = Value
}

/**
  * Kinds of executables that might be generated:
  *  Native: a native platform app or applet
  *  Task: call a task, execute a shell command (usually)
  *  WfFragment: WDL workflow fragment, can included nested if/scatter blocks
  *  WfInputs handle workflow inputs for unlocked workflows
  *  WfOutputs: evaluate workflow outputs
  *  WorkflowOutputReorg: move intermediate result files to a subdirectory.
  */
sealed trait ExecutableKind
case class ExecutableKindNative(executableType: ExecutableType.ExecutableType,
                                id: Option[String] = None,
                                name: Option[String] = None,
                                project: Option[String] = None,
                                path: Option[String] = None)
    extends ExecutableKind
case object ExecutableKindApplet extends ExecutableKind

/**
  * An applet that executes a workflow fragment.
  * @param call name of the call made in the fragment
  * @param blockPath path to the block represented by this fragment
  * @param inputs mapping of input name to type, where names are encoded
  *               such that any dots are replaced with '\_\_\_'
  * @param scatterChunkSize maximum number of scatter jobs that can be
  *                         run at the same time
  */
case class ExecutableKindWfFragment(call: Option[String],
                                    blockPath: Vector[Int],
                                    inputs: Map[DxName, Type],
                                    scatterChunkSize: Option[Int])
    extends ExecutableKind
case class ExecutableKindWfInputs(blockPath: Vector[Int]) extends ExecutableKind
// Output - default and custom reorg
case class ExecutableKindWfOutputs(blockPath: Vector[Int], level: Level.Level)
    extends ExecutableKind
case object ExecutableKindWfCustomReorgOutputs extends ExecutableKind
// Reorg - default and custom reorg
case object ExecutableKindWorkflowOutputReorg extends ExecutableKind
case class ExecutableKindWorkflowCustomReorg(id: String) extends ExecutableKind

object ExecutableKind {
  def getCommand(kind: ExecutableKind): Option[String] = {
    kind match {
      case _: ExecutableKindWfInputs          => Some("Inputs")
      case _: ExecutableKindWfOutputs         => Some("Outputs")
      case ExecutableKindWfCustomReorgOutputs => Some("CustomReorgOutputs")
      case ExecutableKindWorkflowOutputReorg  => Some("OutputReorg")
      case _                                  => None
    }
  }

  def toString(kind: ExecutableKind): String = {
    kind match {
      case _: ExecutableKindNative               => "Native"
      case _: ExecutableKindWfFragment           => "Fragment"
      case ExecutableKindApplet                  => "Task"
      case _: ExecutableKindWfInputs             => "Inputs"
      case _: ExecutableKindWfOutputs            => "Outputs"
      case ExecutableKindWfCustomReorgOutputs    => "Reorg outputs"
      case ExecutableKindWorkflowOutputReorg     => "Output Reorg"
      case ExecutableKindWorkflowCustomReorg(id) => s"Custom reorg ${id}"
    }
  }
}

/**
  * Marker trait for runtime requirements.
  */
trait RuntimeRequirement

/**
  * An app or applet.
  * @param name name of application
  * @param inputs input arguments
  * @param outputs output arguments
  * @param instanceType a platform instance name
  * @param container is docker used? if so, what image
  * @param kind kind of application: task, scatter, ...
  * @param document task definition
  * @param attributes additional applet metadata
  * @param requirements runtime resource requirements
  */
case class Application(name: String,
                       inputs: Vector[Parameter],
                       outputs: Vector[Parameter],
                       instanceType: InstanceType,
                       container: ContainerImage,
                       kind: ExecutableKind,
                       document: SourceCode,
                       attributes: Vector[CallableAttribute] = Vector.empty,
                       requirements: Vector[RuntimeRequirement] = Vector.empty,
                       tags: Set[String] = Set.empty,
                       properties: Map[String, String] = Map.empty,
                       staticFileDependencies: Set[String] = Set.empty)
    extends Callable {
  def inputVars: Vector[Parameter] = inputs
  def outputVars: Vector[Parameter] = outputs
}

/**
  * An input to a stage. Could be empty, a wdl constant,
  * a link to an output variable from another stage,
  * or a workflow input.
  */
sealed trait StageInput
case object StageInputEmpty extends StageInput
case class StageInputStatic(value: Value) extends StageInput
case class StageInputStageLink(stageId: DxWorkflowStage, param: Parameter) extends StageInput
case class StageInputWorkflowLink(param: Parameter) extends StageInput
case class StageInputArray(stageInputs: Vector[StageInput]) extends StageInput
// TODO: the following is not supported directly by DNAnexus, but is needed by at least one
//  CWL conformance test (dynresreq-workflow-stepdefault). We could implement it by making the main
//  stage input optional and using an additional step input to holds the default value.
//case class StageInputWithDefault(stageInput: StageInput, default: Value)

/**
  * A stage can call an application or a workflow.
  * @note the description may contain dots, parentheses, and other special
  * symbols. It is shown to the user on the UI. The dxStage.id is unique
  * across the workflow.
  */
case class Stage(description: String,
                 dxStage: DxWorkflowStage,
                 calleeName: String,
                 inputs: Vector[StageInput],
                 outputs: Vector[Parameter])

/**
  * A workflow output is linked to the stage that generated it.
  *
  * If [level] is SubWorkflow, then a workflow matches part of a
  * WDL workflow, it is not a first class citizen. It is compiled
  * into a hidden dx:workflow.
  */
object Level extends Enum {
  type Level = Value
  val Top, Sub = Value
}

case class Workflow(name: String,
                    inputs: Vector[(Parameter, StageInput)],
                    outputs: Vector[(Parameter, StageInput)],
                    stages: Vector[Stage],
                    document: SourceCode,
                    locked: Boolean,
                    level: Level.Level,
                    attributes: Vector[CallableAttribute] = Vector.empty,
                    tags: Set[String] = Set.empty,
                    properties: Map[String, String] = Map.empty,
                    staticFileDependencies: Set[String] = Set.empty)
    extends Callable {
  def inputVars: Vector[Parameter] = inputs.map(_._1)
  def outputVars: Vector[Parameter] = outputs.map(_._1)
}
