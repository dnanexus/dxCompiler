package dx.core.ir

import dx.core.Constants

import scala.util.matching.Regex

object DxName {
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

  def lookup[T](key: DxName, map: Map[DxName, T]): Option[(DxName, T)] = {
    // the key may be prefixed one or more namespaces - we start with the full name
    // and remove the namespaces successively (left to right) until we find a match
    key
      .dropNamespaceIter()
      .map { key =>
        map.get(key).map(value => (key, value))
      }
      .find {
        case Some(_) => true
        case _       => false
      }
      .flatten
  }
}

/**
  * A name that must conform to the character restrictions for DNAnexus input and output parameter
  * names (`[a-zA-Z0-9_]`).
  *
  * Provides methods for converting the name to its encoded or decoded representation, where
  * "encoded" means that it only contains allowed characters and "decoded" means the language-native
  * form of the name.
  *
  * A name consists of an optional stage prefix, zero or more namespaces, an identifier, and an
  * optional suffix. Name components are delimited by a language-specific delimiter
  * (`namespaceDelim`). When the name is encoded, the namespace  delimiter is encoded as `"___"`,
  * and other disallowed characters are encoded in a language-specific manner. The stage prefix and
  * suffix are added as-is after the encoding/decoding of the name components. We also disallow
  * parts other than the first that begin with an '_' or parts other than the last that end with an
  * underscore, e.g. "foo_._bar".
  *
  * @example {{{
  *   stage prefix: "stage-1"
  *   namespaces: ["foo", "bar"]
  *   identifier: "baz_1"
  *   suffix: "___dxfiles"
  *   encoded: "stage-1.foo___bar___baz_1___dxfiles"
  *   WDL decoded: "stage-1.foo.bar.baz_1___dxfiles"
  *   CWL decoded: "stage-1.foo/bar/baz_1___dxfiles"
  * }}}
  */
abstract class DxName(private var encodedParts: Option[Vector[String]],
                      private var decodedParts: Option[Vector[String]],
                      val stage: Option[String],
                      val suffix: Option[String])
    extends Ordered[DxName] {

  protected def illegalDecodedSequencesRegex: Option[Regex] = None

  assert(
      encodedParts.isEmpty || decodedParts.isEmpty || encodedParts.get.size == decodedParts.get.size,
      s"""encoded and decoded parts are not the same size:
         |${encodedParts.get.size} != ${decodedParts.get.size}""".stripMargin.replaceAll("\n", " ")
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
  decodedParts.foreach { d =>
    assert(d.nonEmpty, s"there must be at least one decoded part")
    d.zipWithIndex.foreach {
      case (part, idx) =>
        assert(part.trim.nonEmpty, s"one or more decoded parts is empty")
        illegalDecodedSequencesRegex.foreach { r =>
          r.findFirstIn(part).map { seq =>
            throw new Exception(s"decoded part ${part} contains illegal sequences '${seq}'")
          }
        }
    }
  }
  stage.foreach { stg =>
    assert(stg.startsWith("stage-"), s"not a valid stage ID ${stg}")
  }

  override def compare(that: DxName): Int = {
    encoded.compare(that.encoded)
  }

  override def hashCode(): Int = {
    encoded.hashCode
  }

  override def equals(obj: Any): Boolean = {
    obj match {
      case that: DxName => compare(that) == 0
      case _            => false
    }
  }

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

  protected def namespaceDelimEncoded: Option[String]

  /**
    * The encoded form of this DxName.
    */
  def encoded: String = {
    val name = getEncodedParts match {
      case Vector(id) => namespaceDelimEncoded.getOrElse("") + id
      case _ if namespaceDelimEncoded.isEmpty =>
        throw new Exception("this name does not allow namespaces")
      case v => namespaceDelimEncoded.get + v.mkString(namespaceDelimEncoded.get)
    }
    (stage, suffix) match {
      case (Some(stg), Some(suf)) => s"${stg}.${name}${suf}"
      case (Some(stg), None)      => s"${stg}.${name}"
      case (None, Some(suf))      => s"${name}${suf}"
      case (None, None)           => name
    }
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

  protected def namespaceDelim: Option[String]

  /**
    * The decoded form of this DxName.
    */
  def decoded: String = {
    val name = getDecodedParts match {
      case Vector(id) => id
      case _ if namespaceDelim.isEmpty =>
        throw new Exception("this name does not allow namespaces")
      case v => v.mkString(namespaceDelim.get)
    }
    (stage, suffix) match {
      case (Some(stg), Some(suf)) => s"${stg}.${name}${suf}"
      case (Some(stg), None)      => s"${stg}.${name}"
      case (None, Some(suf))      => s"${name}${suf}"
      case (None, None)           => name
    }
  }

  def decodedIdentifier: String = {
    getDecodedParts.last
  }

  protected def create(encodedParts: Option[Vector[String]] = None,
                       decodedParts: Option[Vector[String]] = None,
                       stage: Option[String],
                       suffix: Option[String]): DxName

  def withStage(newStage: String): DxName = {
    if (stage.isDefined) {
      throw new Exception(s"${this} already has a stage")
    }
    create(encodedParts, decodedParts, Some(newStage), suffix)
  }

  def dropStage: DxName = {
    create(encodedParts, decodedParts, None, suffix)
  }

  /**
    * Copies this DxName and set suffix to `newSuffix`. Throws an exception
    * if suffix is already defined.
    */
  def withSuffix(newSuffix: String): DxName = {
    if (suffix.isDefined) {
      throw new Exception(s"${this} already has a suffix")
    }
    create(encodedParts, decodedParts, stage, Some(newSuffix))
  }

  /**
    * Copies this DxName and adds `newSuffix` after the existing suffix, if any.
    */
  def addSuffix(suffixToAdd: String): DxName = {
    create(encodedParts,
           decodedParts,
           stage,
           suffix.map(suf => s"${suf}${suffixToAdd}").orElse(Some(suffixToAdd)))
  }

  /**
    * Copies this DxName with `suffixToDrop` removed from the existing
    * suffix (if any). If `suffixToDrop` is `None` or equal to the
    * existing suffix, then the existing suffix is set to `None`.
    */
  def dropSuffix(suffixToDrop: Option[String] = None): DxName = {
    val newSuffix = (suffix, suffixToDrop) match {
      case (Some(s1), Some(s2)) if s1 == s2        => None
      case (Some(s1), Some(s2)) if s1.endsWith(s2) => Some(s1.dropRight(s2.length))
      case (s1, Some(s2)) =>
        throw new Exception(s"existing suffix ${s1} does not end with ${s2}")
      case _ => None
    }
    create(encodedParts, decodedParts, stage, newSuffix)
  }

  def endsWith(dxName: DxName): Boolean = {
    decoded.endsWith(dxName.decoded)
  }

  /**
    * Creates a new DxName with `ns` inserted as the first namespace component.
    */
  def pushDecodedNamespace(ns: String): DxName = {
    create(decodedParts = Some(Vector(ns) ++ getDecodedParts), stage = stage, suffix = suffix)
  }

  def pushDecodedNamespaces(namespaces: Vector[String]): DxName = {
    create(decodedParts = Some(namespaces ++ getDecodedParts), stage = stage, suffix = suffix)
  }

  /**
    * Removes the first namespace part and returns it as a String along with the new name.
    * Throws an Exception if there are no namespace parts.
    */
  def popDecodedNamespace(keepStage: Boolean = false): (String, DxName) = {
    val decoded = getDecodedParts
    if (decoded.size <= 1) {
      throw new Exception("cannot pop identifier unless there is more than one part")
    }
    val newName = create(decodedParts = Some(decoded.drop(1)),
                         stage = Option.when(keepStage)(stage).flatten,
                         suffix = suffix)
    (decoded.head, newName)
  }

  def dropNamespaceIter(keepStage: Boolean = false): Iterator[DxName] = {
    Iterator.unfold[DxName, Option[DxName]](Some(this)) {
      case Some(n) if n.numParts > 1 => Some((n, Some(n.popDecodedNamespace(keepStage)._2)))
      case Some(n)                   => Some((n, None))
      case _                         => None
    }
  }

  /**
    * Returns a new DxName with only the identifier part and the suffix.
    */
  def dropNamespaces(keepStage: Boolean = false): DxName = {
    if (numParts > 1) {
      create(
          encodedParts = encodedParts.map(v => Vector(v.last)),
          decodedParts = decodedParts.map(v => Vector(v.last)),
          stage = Option.when(keepStage)(stage).flatten,
          suffix = suffix
      )
    } else {
      this
    }
  }

  /**
    * Removes the identifier (the last part) and returns it as a String along with the new name.
    * Throws an Exception if there are no namespace parts. The suffix is retained in the new name
    * only if `keepSuffix=true`.
    */
  def popDecodedIdentifier(keepSuffix: Boolean = false): (DxName, String) = {
    val decoded = getDecodedParts
    if (decoded.size <= 1) {
      throw new Exception("cannot pop identifier unless there is more than one part")
    }
    val newName = create(decodedParts = Some(decoded.dropRight(1)),
                         stage = stage,
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
  private val stageRegex = "^(?:((?:stage-[^.]+\\.)*stage-[^.]+)\\.)?(.+)$".r

  // the default Regex.split method does not return parts with empty strings -
  // we need those for validation
  implicit class RegexExtensions(regex: Regex) {
    def split(toSplit: CharSequence, limit: Int): Array[String] = {
      regex.pattern.split(toSplit, limit)
    }
  }

  def split(name: String,
            sepRegex: Option[Regex] = None): (Vector[String], Option[String], Option[String]) = {
    def splitSuffix(s: String): (String, Option[String]) = {
      suffixes
        .collectFirst {
          case suffix if s.endsWith(suffix) =>
            val (prefix, firstSuffix) = splitSuffix(s.dropRight(suffix.length))
            (prefix, firstSuffix.map(f => s"${f}${suffix}").orElse(Some(suffix)))
        }
        .getOrElse(s, None)
    }
    def splitParts(s: String): Vector[String] = {
      sepRegex.map(_.split(s, Int.MaxValue).toVector).getOrElse(Vector(s))
    }
    val (stage, rest) = name match {
      case stageRegex(stage, rest) => (Option(stage), rest)
    }
    val (prefix, suffix) = splitSuffix(rest)
    val parts = splitParts(prefix)
    val validparts = parts.zipWithIndex.foldLeft(Vector.empty[String]) {
      case (accu, (part, idx)) if !part.isEmpty || (idx != 0) => //&& idx != parts.size - 1) =>
        accu ++ Vector(part)
      case (accu, (_, _)) => accu
    }
    (validparts, stage, suffix)
  }
}

object SimpleDxName extends DxNameFactory {
  override def fromEncodedName(name: String): SimpleDxName = {
    val (parts, stage, suffix) = DxNameFactory.split(name)
    new SimpleDxName(parts, stage, suffix)
  }

  override def fromDecodedName(name: String): SimpleDxName = {
    val (parts, stage, suffix) = DxNameFactory.split(name)
    new SimpleDxName(parts, stage, suffix)
  }

  def fromSourceName(identifier: String, suffix: Option[String] = None): SimpleDxName = {
    new SimpleDxName(Vector(identifier), stage = None, suffix = suffix)
  }
}

/**
  * A DxName that does not perform encoding or decoding - the encoded and
  * decoded names must be equal.
  */
class SimpleDxName(parts: Vector[String],
                   stage: Option[String] = None,
                   suffix: Option[String] = None,
                   override val namespaceDelim: Option[String] = None)
    extends DxName(Some(parts), Some(parts), stage, suffix) {
  override protected def create(encodedParts: Option[Vector[String]],
                                decodedParts: Option[Vector[String]],
                                stage: Option[String],
                                suffix: Option[String]): DxName = {
    (encodedParts, decodedParts) match {
      case (Some(e), Some(d)) if DxName.compareStringVectors(e, d) != 0 =>
        throw new Exception(s"encoded and decoded parts must be the same: ${e} != ${d}")
      case _ =>
        new SimpleDxName(encodedParts.orElse(decodedParts).get, stage, suffix, namespaceDelim)
    }
  }
}
