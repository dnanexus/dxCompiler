package dx.core.ir

import dx.api.DxWorkflowStage
import dx.core.Constants
import dx.core.ir.RunSpec.{ContainerImage, InstanceType}
import dx.util.Enum

import scala.util.matching.Regex

object DxName {
  val NamespaceDelimEncoded = "___"
  // characters that need to be encoded
  val disallowedCharsRegex: Regex = "([^a-zA-Z0-9_])".r

  def compareStringVectors(a: Vector[String], b: Vector[String]): Int = {
    a.zip(b)
      .iterator
      .map { case (i, j) => i.compare(j) }
      .collectFirst {
        case x if x != 0 => x
      }
      .getOrElse(a.size.compare(b.size))
  }
}

/**
  * A name that must conform to the character restrictions for DNAnexus input
  * and output parameter names (`[a-zA-Z0-9_]`).
  *
  * Provides methods for converting the name to its encoded or decoded
  * representation, where "encoded" means that it only contains allowed characters
  * and "decoded" means the language-native form of the name.
  *
  * A name consists of zero or more namespaces, followed by an identifier, followed
  * by an optional suffix. Name components are delimited by a language-specific
  * delimiter (`namespaceDelim`). When the name is encoded, the namespace
  * delimiter is encoded as `"___"`, and other disallowed characters are encoded in
  * a language-specific manner. The suffix is added as-is after the encoding/decoding
  * of the prefix, so it must not contain any disallowed characters.  We also disallow
  * parts other than the first that begin with an '_' or parts other than the last that
  * end with an underscore, e.g. "foo_._bar".
  *
  * @example {{{
  *   namespaces: ["foo", "bar"]
  *   identifier: "baz_1"
  *   suffix: "___dxfiles"
  *   encoded: "foo___bar___baz_1___dxfiles"
  *   WDL decoded: "foo.bar.baz_1___dxfiles"
  *   CWL decoded: "foo/bar/baz_1___dxfiles"
  * }}}
  */
abstract class DxName(private var encodedParts: Option[Vector[String]],
                      private var decodedParts: Option[Vector[String]],
                      val suffix: Option[String])
    extends Ordered[DxName] {
  assert(
      encodedParts.isEmpty || decodedParts.isEmpty || encodedParts.get.size == decodedParts.get.size,
      s"encoded and decoded parts are not the same size: ${encodedParts.get.size} != ${decodedParts.get.size}"
  )
  encodedParts.foreach { e =>
    assert(e.nonEmpty, s"there must be at least one encoded part")
    e.foreach { part =>
      assert(part.trim.nonEmpty, s"one or more encoded parts is empty")
      DxName.disallowedCharsRegex.findFirstIn(part).map { c =>
        throw new Exception(s"encoded part ${part} contains illegal character '${c}'")
      }
    }
  }
  protected def illegalDecodedSequencesRegex: Option[Regex] = None
  decodedParts.foreach { d =>
    assert(d.nonEmpty, s"there must be at least one decoded part")
    d.zipWithIndex.foreach {
      case (part, idx) =>
        assert(part.trim.nonEmpty, s"one or more decoded parts is empty")
        assert(idx == 0 || !part.startsWith("_"), "only the first part may start with '_'")
        assert(
            !part.endsWith("_") || (idx == d.size - 1 && !suffix.exists(_.startsWith("_"))),
            "only the last part may end with '_', and only if the suffix does not start with '_'"
        )
        illegalDecodedSequencesRegex.foreach { r =>
          r.findFirstIn(part).map { seq =>
            throw new Exception(s"decoded part ${part} contains illegal sequences ${seq}")
          }
        }
    }
  }

  protected def namespaceDelim: Option[String]

  override def compare(that: DxName): Int = {
    val cmp = if (encodedParts.isDefined && that.encodedParts.isDefined) {
      DxName.compareStringVectors(encodedParts.get, that.encodedParts.get)
    } else if (decodedParts.isDefined && that.decodedParts.isDefined) {
      DxName.compareStringVectors(decodedParts.get, that.decodedParts.get)
    } else {
      decoded.compare(that.decoded)
    }
    if (cmp != 0) {
      cmp
    } else {
      (suffix, that.suffix) match {
        case (None, None)         => 0
        case (Some(s1), Some(s2)) => s1.compare(s2)
        case (Some(_), _)         => 1
        case _                    => -1
      }
    }
  }

  override def hashCode(): Int = {
    encoded.hashCode
  }

  // subclasses must override equals - comparison with
  // another name of a different subclass must return false
  override def equals(obj: Any): Boolean = ???

  def numParts: Int = {
    encodedParts.map(_.size).orElse(decodedParts.map(_.size)).get
  }

  protected def encodePart(part: String): String = part

  def getEncodedParts: Vector[String] = {
    encodedParts
      .getOrElse {
        val encoded = decodedParts.get.map(encodePart)
        encodedParts = Some(encoded)
        encoded
      }
  }

  /**
    * The encoded form of this DxName.
    */
  def encoded: String = {
    val prefix = getEncodedParts.mkString(DxName.NamespaceDelimEncoded)
    suffix.map(suffix => s"${prefix}${suffix}").getOrElse(prefix)
  }

  protected def decodePart(part: String): String = part

  def getDecodedParts: Vector[String] = {
    decodedParts
      .getOrElse {
        val decoded = encodedParts.get.map(decodePart)
        decodedParts = Some(decoded)
        decoded
      }
  }

  /**
    * The decoded form of this DxName.
    */
  def decoded: String = {
    val prefix = getDecodedParts match {
      case Vector(id) => id
      case _ if namespaceDelim.isEmpty =>
        throw new Exception("this name does not allow namespaces")
      case v => v.mkString(namespaceDelim.get)
    }
    suffix.map(suffix => s"${prefix}${suffix}").getOrElse(prefix)
  }

  protected def create(encodedParts: Option[Vector[String]] = None,
                       decodedParts: Option[Vector[String]] = None,
                       suffix: Option[String] = None): DxName

  /**
    * Copies this DxName and set suffix to `newSuffix`. Throws an exception
    * if suffix is already defined.
    */
  def withSuffix(newSuffix: String): DxName = {
    if (suffix.isDefined) {
      throw new Exception(s"${this} already has a suffix")
    }
    create(encodedParts, decodedParts, Some(newSuffix))
  }

  /**
    * Copies this DxName and sets suffix to `None`. Throws an exception
    * if suffix is already `None`.
    */
  def dropSuffix: DxName = {
    if (suffix.isEmpty) {
      throw new Exception(s"${this} does not have a suffix")
    }
    create(encodedParts, decodedParts, None)
  }

  def endsWith(dxName: DxName): Boolean = {
    decoded.endsWith(dxName.decoded)
  }

  /**
    * Creates a new DxName with `ns` inserted as the first namespace component.
    */
  def pushDecodedNamespace(ns: String): DxName = {
    create(decodedParts = Some(Vector(ns) ++ getDecodedParts), suffix = suffix)
  }

  def popDecodedIdentifier(keepSuffix: Boolean = false): (DxName, String) = {
    val decoded = getDecodedParts
    if (decoded.size <= 1) {
      throw new Exception("cannot pop identifier unless there is more than one part")
    }
    val newName = create(decodedParts = Some(decoded.dropRight(1)),
                         suffix = Option.when(keepSuffix)(suffix).flatten)
    (newName, decoded.last)
  }

  override def toString: String = decoded
}

trait DxNameFactory {
  def fromEncodedName(name: String): DxName

  def fromDecodedName(name: String): DxName
}

object DxNameFactory {
  // standard suffixes to parse
  val suffixes = Set(Constants.ComplexValueKey, Constants.FlatFilesSuffix)

  // the default Regex.split method does not return parts with empty strings -
  // we need those for validation
  implicit class RegexExtensions(regex: Regex) {
    def split(toSplit: CharSequence, limit: Int): Array[String] = {
      regex.pattern.split(toSplit, limit)
    }
  }

  def split(name: String, sepRegex: Option[Regex] = None): (Vector[String], Option[String]) = {
    def splitParts(s: String): Vector[String] = {
      sepRegex.map(_.split(s, Int.MaxValue).toVector).getOrElse(Vector(s))
    }
    suffixes
      .collectFirst {
        case suffix if name.endsWith(suffix) =>
          (splitParts(name.dropRight(suffix.length)), Some(suffix))
      }
      .getOrElse((splitParts(name), None))
  }
}

object SimpleDxName extends DxNameFactory {
  override def fromEncodedName(name: String): SimpleDxName = {
    val (parts, suffix) = DxNameFactory.split(name)
    new SimpleDxName(parts, suffix)
  }

  override def fromDecodedName(name: String): SimpleDxName = {
    val (parts, suffix) = DxNameFactory.split(name)
    new SimpleDxName(parts, suffix)
  }

  def fromSourceName(identifier: String, suffix: Option[String] = None): SimpleDxName = {
    new SimpleDxName(Vector(identifier), suffix)
  }
}

/**
  * A DxName that does not perform encoding or decoding - the encoded and
  * decoded names must be equal.
  */
class SimpleDxName(parts: Vector[String],
                   suffix: Option[String] = None,
                   override val namespaceDelim: Option[String] = None)
    extends DxName(Some(parts), Some(parts), suffix) {

  override def equals(obj: Any): Boolean = {
    obj match {
      case that: SimpleDxName => compare(that) == 0
      case _                  => false
    }
  }

  override protected def create(encodedParts: Option[Vector[String]],
                                decodedParts: Option[Vector[String]],
                                suffix: Option[String]): DxName = {
    (encodedParts, decodedParts) match {
      case (Some(e), Some(d)) if DxName.compareStringVectors(e, d) != 0 =>
        throw new Exception(s"encoded and decoded parts must be the same: ${e} != ${d}")
      case _ =>
        new SimpleDxName(encodedParts.orElse(decodedParts).get, namespaceDelim, suffix)
    }
  }
}

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
