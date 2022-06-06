package dx.core.languages.wdl

import dx.core.ir.{DxName, DxNameFactory}

import scala.util.matching.Regex

object WdlDxName extends DxNameFactory {
  val NamespaceDelim = "."
  private val NamespaceDelimRegex = "\\.".r
  private val NamespaceDelimEncodedRegex: Regex = DxName.NamespaceDelimEncoded.r

  override def fromEncodedName(name: String): WdlDxName = {
    val (parts, stage, suffix) = DxNameFactory.split(name, Some(NamespaceDelimEncodedRegex))
    new WdlDxName(encodedParts = Some(parts), stage = stage, suffix = suffix)
  }

  def empty: WdlDxName = new WdlDxName()

  override def fromDecodedName(name: String): WdlDxName = {
    val (parts, stage, suffix) = DxNameFactory.split(name, Some(NamespaceDelimRegex))
    new WdlDxName(decodedParts = Some(parts), stage = stage, suffix = suffix)
  }

  /**
    * Creates a WdlDxName from an identifier that does not contain any disallowed characters.
    */
  def fromSourceName(identifier: String,
                     namespace: Option[String] = None,
                     suffix: Option[String] = None): WdlDxName = {
    new WdlDxName(decodedParts = Some(namespace.toVector ++ Vector(identifier)),
                  stage = None,
                  suffix = suffix)
  }
}

/**
  * Simplified decoding for parameter names that conform to DNAnexus character
  * restrictions, with the possible exception of containg dots.
  */
class WdlDxName(encodedParts: Option[Vector[String]] = None,
                decodedParts: Option[Vector[String]] = None,
                stage: Option[String] = None,
                suffix: Option[String] = None)
    extends DxName(encodedParts, decodedParts, stage, suffix) {

  override protected def illegalDecodedSequencesRegex: Option[Regex] =
    Some(DxName.disallowedCharsRegex)

  override protected val namespaceDelim: Option[String] = Some(WdlDxName.NamespaceDelim)

  override protected def create(encodedParts: Option[Vector[String]] = None,
                                decodedParts: Option[Vector[String]] = None,
                                stage: Option[String],
                                suffix: Option[String]): WdlDxName = {
    new WdlDxName(encodedParts, decodedParts, stage, suffix)
  }
}
