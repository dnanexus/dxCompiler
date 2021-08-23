package dx.core.ir

import dx.api.DxWorkflowStage
import dx.core.ir.RunSpec.{InstanceType, ContainerImage}
import dx.util.Enum

trait ParameterAttribute

/**
  * Compile time representation of a variable. Used also as an applet argument.
  *
  * The fullyQualifiedName could contains dots. However dx does not allow dots
  * in applet/workflow arugment names, this requires some kind of transform.
  *
  * The attributes are used to encode DNAx applet input/output specification
  * fields, such as {help, suggestions, patterns}.
  *
  * @param name parameter name
  * @param dxType parameter data type
  * @param defaultValue default value
  * @param attributes metadata
  */
case class Parameter(
    name: String,
    dxType: Type,
    defaultValue: Option[Value] = None,
    attributes: Vector[ParameterAttribute] = Vector.empty
) {
  // dx does not allow dots, slashes, etc in variable names, so we encode them.
  // TODO: check for collisions that are created this way.
  def dxName: String = Parameter.encodeName(name)
}

object Parameter {
  val DefaultDelim = "."
  val ComplexValueKey = "___"

  // character sequences that may not appear in a non-encoded name
  private val illegalSeqsRegexp = "__|\\s+".r
  // characters that need to be encoded
  private val encodeCharsRegexp = "([^a-zA-Z0-9_])".r
  // character sequences that need to be decoded
  private val decodeSeqsRegexp = "(___)|__(\\d+)__".r

  private def makeNewName(name: String,
                          starts: Vector[Int],
                          ends: Vector[Int],
                          delims: Vector[String]): String = {
    (Vector(0) ++ starts)
      .zip(ends ++ Vector(name.length))
      .map {
        case (start, end) if start == end => ""
        case (start, end)
            if (start > 0 && name.charAt(start) == '_')
              || (end < name.length && name.charAt(end - 1) == '_') =>
          println(start, end, name.substring(start, end))
          throw new Exception(
              s"illegal name ${name}: '_' must not be adjacent to a non-alphanumeric character"
          )
        case (start, end) => name.substring(start, end)
      }
      .zipAll(delims, "", "")
      .map {
        case (part, delim) => s"${part}${delim}"
      }
      .mkString
  }

  /**
    * Converts disallowed characters in parameter names to underscores.
    * DNAnexus only allows [a-zA-Z0-9_] in parameter names. We also
    * disallow:
    *   * whitespace - "foo bar"
    *   * consecutive underscores - "foo\_\_bar"
    *   * name segments that begin or end in underscore - "foo_._bar"
    * @param name parameter name
    * @return the encoded parameter name, with "." replaced with "\_\_\_"
    *         and any other illegal character replaced with "\_\_ord(x)\_\_",
    *         where `ord(x)` is the ordinal value of the illegal character.
    */
  def encodeName(name: String): String = {
    if (name.trim().isEmpty) {
      throw new Exception("empty name")
    }
    if (name.endsWith(".")) {
      throw new Exception(s"parameter name ${name} ends with a '.'")
    }
    // some special parameter names end with '___' - add this back on after decoding
    val (encodeName, endsWithComplexValueKey) = if (name.endsWith(ComplexValueKey)) {
      (name.dropRight(ComplexValueKey.length), true)
    } else {
      (name, false)
    }
    illegalSeqsRegexp.findFirstIn(encodeName).foreach { c =>
      throw new Exception(s"parameter name contains illegal character sequence '${c}'")
    }
    val encoded = encodeCharsRegexp.findAllMatchIn(encodeName).toVector match {
      case Vector() => name
      case matches =>
        val (starts, ends, delims) = matches.map { m =>
          val newDelim = m.group(1) match {
            case DefaultDelim => ComplexValueKey
            case c if c.length() == 1 =>
              s"__${c.charAt(0).toInt}__"
            case other =>
              throw new Exception(s"unexpected multi-character delimiter ${other}")
          }
          (m.end, m.start, newDelim)
        }.unzip3
        makeNewName(encodeName, starts, ends, delims)
    }
    if (endsWithComplexValueKey) {
      s"${encoded}${ComplexValueKey}"
    } else {
      encoded
    }
  }

  /**
    * Simplified encoding for parameter  names that conform to DNAnexus character
    * restrictions, with the possible exception of containg dots.
    * @param name parameter name
    * @return encoded parameter name
    */
  def encodeDots(name: String): String = {
    if (name.trim().isEmpty) {
      throw new Exception("empty name")
    }
    if (name.endsWith(".")) {
      throw new Exception(s"parameter name ${name} ends with a '.'")
    }
    name.split("\\.").toVector match {
      case Vector(_) => name
      case parts     =>
        // check that none of the individual words start/end with '_',
        // which will cause problems when decoding
        if (parts.dropRight(1).exists(_.endsWith("_")) || parts.drop(1).exists(_.startsWith("_"))) {
          throw new Exception(s"cannot encode value ${name} - cannot have '_' next to '.'")
        }
        parts.mkString(ComplexValueKey)
    }
  }

  /**
    * Decodes a name that was previously encoded using `encodeName`.
    * @param name encoded parameter name
    * @return parameter name
    */
  def decodeName(name: String): String = {
    // some special parameter names end with '___' - add this back on after decoding
    val (decodeName, endsWithComplexValueKey) = if (name.endsWith(ComplexValueKey)) {
      (name.dropRight(ComplexValueKey.length), true)
    } else {
      (name, false)
    }
    if (decodeName.trim().isEmpty) {
      throw new Exception("empty name")
    }
    encodeCharsRegexp.findFirstIn(decodeName).foreach { c =>
      throw new Exception(s"encoded value ${name} contains illegal character '${c}''")
    }
    decodeSeqsRegexp.findAllMatchIn(decodeName).toVector match {
      case Vector() => name
      case matches =>
        val (starts, ends, delims) = matches.map { m =>
          val newDelim = (m.group(1), m.group(2)) match {
            case (ComplexValueKey, null) => DefaultDelim
            case (null, ord)             => ord.toInt.toChar.toString
            case other =>
              throw new Exception(s"unexpected delimiter ${other}")
          }
          (m.end, m.start, newDelim)
        }.unzip3
        val newName = makeNewName(decodeName, starts, ends, delims)
        if (endsWithComplexValueKey) {
          s"${newName}${ComplexValueKey}"
        } else {
          newName
        }
    }
  }

  /**
    * Simplified decoding for parameter names that conform to DNAnexus character
    * restrictions, with the possible exception of containg dots.
    * @param name encoded parameter name
    * @return parameter name
    */
  def decodeDots(name: String): String = {
    if (name.trim().isEmpty) {
      throw new Exception("empty name")
    }
    if (name.contains(".")) {
      throw new Exception(s"encoded name ${name} contains '.'")
    }
    // some special parameter names end with '___' - add this back on after decoding
    val (decodeName, endsWithComplexValueKey) = if (name.endsWith(ComplexValueKey)) {
      (name.dropRight(ComplexValueKey.length), true)
    } else {
      (name, false)
    }
    decodeName.split(ComplexValueKey).toVector match {
      case Vector(_) => name
      case parts     =>
        // check that none of the individual words start/end with '_'
        if (parts.dropRight(1).exists(_.endsWith("_")) || parts.drop(1).exists(_.startsWith("_"))) {
          throw new Exception(
              s"cannot decode value ${name} - more than three consecutive underscores"
          )
        }
        val decoded = parts.mkString(DefaultDelim)
        if (endsWithComplexValueKey) {
          // some special parameter names end with '___' - add this back on after decoding
          s"${decoded}${ComplexValueKey}"
        } else {
          decoded
        }
    }
  }
}

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
  * @param call names of calls made in the fragment
  * @param blockPath path to the block represented by this fragment
  * @param inputs mapping of input name to type, where names are encoded
  *               such that any dots are replaced with '\_\_\_'
  * @param scatterChunkSize maximum number of scatter jobs that can be
  *                         run at the same time
  */
case class ExecutableKindWfFragment(call: Option[String],
                                    blockPath: Vector[Int],
                                    inputs: Map[String, Type],
                                    scatterChunkSize: Option[Int])
    extends ExecutableKind
case class ExecutableKindWfInputs(blockPath: Vector[Int]) extends ExecutableKind
// Output - default and custom reorg
case class ExecutableKindWfOutputs(blockPath: Vector[Int]) extends ExecutableKind
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
case object EmptyInput extends StageInput
case class StaticInput(value: Value) extends StageInput
case class LinkInput(stageId: DxWorkflowStage, paramName: String) extends StageInput
case class WorkflowInput(param: Parameter) extends StageInput
case class ArrayInput(stageInputs: Vector[StageInput]) extends StageInput

// A stage can call an application or a workflow.
//
// Note: the description may contain dots, parentheses, and other special
// symbols. It is shown to the user on the UI. The dxStage.id is unique
// across the workflow.
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
  def inputVars: Vector[Parameter] = inputs.map { case (cVar, _)   => cVar }
  def outputVars: Vector[Parameter] = outputs.map { case (cVar, _) => cVar }
}
