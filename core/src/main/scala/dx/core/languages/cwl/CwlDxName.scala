package dx.core.languages.cwl

import dx.core.ir.{DxName, DxNameFactory}

import scala.util.matching.Regex

object CwlDxName extends DxNameFactory {
  val NamespaceDelim = "/"
  private val NamespaceDelimRegex = NamespaceDelim.r
  // Regex for finding exactly three consecutive underscores. When
  // encoding consecutive disallowed characters, we will end up with
  // something like "__42____67__", which we don't want to match as
  // a namespace delimiter.
  private val encodedNamespaceDelimRegex = "((?<=[^_])|^)(___)((?=[^_])|$)".r
  // character sequences that may not appear in a non-encoded name
  private val illegalDecodedSequencesRegex = "/|__|\\s+".r

  private val disallowedDxNameFormat: String = "^[^a-zA-Z_].*"

  override def fromEncodedName(name: String): CwlDxName = {
    val (parts, stage, suffix) = DxNameFactory.split(name, Some(encodedNamespaceDelimRegex))
    new CwlDxName(encodedParts = Some(parts), stage = stage, suffix = suffix)
  }

  override def fromDecodedName(name: String): CwlDxName = {
    val (parts, stage, suffix) = DxNameFactory.split(name, Some(NamespaceDelimRegex))
    new CwlDxName(decodedParts = Some(parts), stage = stage, suffix = suffix)
  }

  /**
    * Creates a ComplexDxName from a parameter name coming directly from a
    * language AST, with optional namespace.
    */
  def fromSourceName(identifier: String,
                     namespace: Option[String] = None,
                     suffix: Option[String] = None): CwlDxName = {
    new CwlDxName(decodedParts = Some(namespace.toVector ++ Vector(identifier)),
                  stage = None,
                  suffix = suffix)
  }
}

class CwlDxName(encodedParts: Option[Vector[String]] = None,
                decodedParts: Option[Vector[String]] = None,
                stage: Option[String] = None,
                suffix: Option[String] = None)
    extends DxName(encodedParts, decodedParts, stage, suffix) {

  override protected def illegalDecodedSequencesRegex: Option[Regex] =
    Some(CwlDxName.illegalDecodedSequencesRegex)

  override protected def namespaceDelim: Option[String] = Some(CwlDxName.NamespaceDelim)

  override protected def create(encodedParts: Option[Vector[String]] = None,
                                decodedParts: Option[Vector[String]] = None,
                                stage: Option[String],
                                suffix: Option[String]): CwlDxName = {
    new CwlDxName(encodedParts, decodedParts, stage, suffix)
  }

  /**
    * The encoded form of CwlDxName.
    */
  override def encoded: String = {
    val name = getEncodedParts.mkString(DxName.NamespaceDelimEncoded) match {
      case n if n.matches(CwlDxName.disallowedDxNameFormat) =>
        DxName.NamespaceDelimEncoded + n
      case n => n
    }
    (stage, suffix) match {
      case (Some(stg), Some(suf)) => s"${stg}.${name}${suf}"
      case (Some(stg), None)      => s"${stg}.${name}"
      case (None, Some(suf))      => s"${name}${suf}"
      case (None, None)           => name
    }
  }
}
