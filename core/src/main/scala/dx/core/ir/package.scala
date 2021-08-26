package dx.core.ir

import dx.api.DxWorkflowStage
import dx.core.Constants
import dx.core.ir.RunSpec.{ContainerImage, InstanceType}
import dx.util.Enum

trait ParameterAttribute

object DxName {
  val Dot = "."
  val DotEncoded = "___"
  // standard suffixes to parse
  val suffixes = Set(Constants.ComplexValueKey, Constants.FlatFilesSuffix)

  def split(name: String): (String, Option[String]) = {
    DxName.suffixes
      .collectFirst {
        case suffix if name.endsWith(suffix) => (name.dropRight(suffix.length), Some(suffix))
      }
      .getOrElse((name, None))
  }
}

/**
  * A name that must conform to the character restrictions for DNAnexus input
  * and output parameter names. May be created from an encoded prefix and/or
  * a decoded prefix, and an optional suffix. The encoded prefix is the
  * DNAnexus-compatible version of the name, and the decoded prefix is the
  * native version of the name. The suffix is added to the end of either prefix
  * without any encoding.
  * @example
  * decodedPrefix = "foo.bar_baz"
  * encodedPrefix = "foo\_\_\_bar_baz"
  * suffix = "\_\_\_dxfiles"
  * encoded = "foo\_\_\_bar_baz\_\_\_dxfiles"
  * decoded = foo.bar_baz\_\_\_dxfiles
  *
  * TODO: write tests for equality and hash lookup
  */
abstract class DxName(private var encodedPrefix: Option[String],
                      private var decodedPrefix: Option[String],
                      val suffix: Option[String])
    extends Ordered[DxName] {
  assert(encodedPrefix.exists(_.trim().nonEmpty) || decodedPrefix.exists(_.trim().nonEmpty),
         "at least one of encodedPrefix, decodedPrefix is required")

  override def compare(that: DxName): Int = {
    val cmp = if (encodedPrefix.isDefined && that.encodedPrefix.isDefined) {
      encodedPrefix.get.compare(that.encodedPrefix.get)
    } else if (decodedPrefix.isDefined && that.decodedPrefix.isDefined) {
      decodedPrefix.get.compare(that.decodedPrefix.get)
    } else {
      decoded.compare(that.decoded)
    }
    if (cmp != 0) {
      cmp
    } else {
      (suffix, that.suffix) match {
        case (None, None)         => 0
        case (Some(s1), Some(s2)) => s1.compare(s2)
        case (_, Some(_))         => 1
        case _                    => -1
      }
    }
  }

  override def hashCode(): Int = {
    decoded.hashCode
  }

  override def equals(obj: Any): Boolean = {
    obj match {
      case that: DxName => compareTo(that) == 0
      case _            => false
    }
  }

  protected def encodePrefix(prefix: String): String

  def encoded: String = {
    val encoded = encodedPrefix.getOrElse {
      val encoded = encodePrefix(decodedPrefix.get)
      encodedPrefix = Some(encoded)
      encoded
    }
    suffix.map(suf => s"${encoded}${suf}").getOrElse(encoded)
  }

  protected def decodePrefix(prefix: String): String

  def decoded: String = {
    val decoded = decodedPrefix.getOrElse {
      val decoded = decodePrefix(encodedPrefix.get)
      decodedPrefix = Some(decoded)
      decoded
    }
    suffix.map(suf => s"${decoded}${suf}").getOrElse(decoded)
  }

  protected def create(encodedPrefix: Option[String] = None,
                       decodedPrefix: Option[String] = None,
                       suffix: Option[String] = None): DxName

  /**
    * Copies this DxName and set suffix to `newSuffix`. Throws an exception
    * if suffix is already defined.
    */
  def withSuffix(newSuffix: String): DxName = {
    if (suffix.isDefined) {
      throw new Exception(s"${this} already has a suffix")
    }
    create(encodedPrefix, decodedPrefix, Some(newSuffix))
  }

  /**
    * Copies this DxName and sets suffix to `None`. Throws an exception
    * if suffix is already `None`.
    */
  def dropSuffix: DxName = {
    if (suffix.isEmpty) {
      throw new Exception(s"${this} does not have a suffix")
    }
    create(encodedPrefix, decodedPrefix, None)
  }

  def endsWith(dxName: DxName): Boolean = {
    decoded.endsWith(dxName.decoded)
  }

  protected def namespaceDelim: String = DxName.Dot

  /**
    * Creates a new DxName with `prefix` added to the current prefix, delimited
    * by `delim`.
    */
  def addDecodedNamespace(prefix: String): DxName = {
    val newDecodedPrefix =
      s"${prefix}${namespaceDelim}${decodedPrefix.getOrElse(decodePrefix(encodedPrefix.get))}"
    create(decodedPrefix = Some(newDecodedPrefix), suffix = suffix)
  }

  def getDecodedNamespace(dxName: DxName): String = {
    val fqn = decoded
    val suffix = s"${namespaceDelim}${dxName.decoded}"
    if (fqn.endsWith(suffix)) {
      fqn.dropRight(fqn.length)
    } else {
      throw new Exception(s"${fqn} does not end with ${suffix}")
    }
  }

  // TODO: for now throw an exception because any usage of this might be wrong -
  //  eventually change this to return `decoded`
  override def toString: String = {
    throw new Exception("DxName.toString")
  }
}

/**
  * A DxName that does not perform encoding or decoding - the encoded and
  * decoded names must be equal.
  */
class SimpleDxName(prefix: String, suffix: Option[String] = None)
    extends DxName(Some(prefix), Some(prefix), suffix) {
  override protected def encodePrefix(prefix: String): String = {
    throw new Exception("unreachable")
  }

  override protected def decodePrefix(prefix: String): String = {
    throw new Exception("unreachable")
  }

  override protected def create(encodedPrefix: Option[String],
                                decodedPrefix: Option[String],
                                suffix: Option[String]): DxName = {
    (encodedPrefix, decodedPrefix) match {
      case (Some(p1), Some(p2)) if p1 == p2 => new SimpleDxName(p1, suffix)
      case _ =>
        throw new Exception(
            s"encoded and decoded names are not equal: ${encodedPrefix} != ${decodedPrefix}"
        )
    }
  }
}

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
  * @param call names of calls made in the fragment
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
case class LinkInput(stageId: DxWorkflowStage, paramName: DxName) extends StageInput
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
