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
  // character sequences that need to be decoded
  private val decodeSequencesRegex = "__(\\d+)__".r

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

  private def makeNewPart(name: String,
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
    * Encodes disallowed characters in a part. In addition to DxName restrictions,
    * we also disallow parts with consecutive underscores, e.g. {{{"foo__bar"}}}.
    * @param part the decoded prefix
    * @return the encoded parameter name with each illegal character
    *         replaced with {{{"__ord(x)__"}}}, where `ord(x)` is the ordinal
    *         value of the illegal character.
    */
  override protected def encodePart(part: String): String = {
    DxName.disallowedCharsRegex.findAllMatchIn(part).toVector match {
      case Vector() => part
      case matches =>
        val (starts, ends, delims) = matches.map { m =>
          val newDelim = m.group(1) match {
            case c if c.length() == 1 => s"__${c.charAt(0).toInt}__"
            case other =>
              throw new Exception(s"unexpected multi-character delimiter ${other}")
          }
          (m.end, m.start, newDelim)
        }.unzip3
        makeNewPart(part, starts, ends, delims)
    }
  }

  override protected def decodePart(part: String): String = {
    CwlDxName.decodeSequencesRegex.findAllMatchIn(part).toVector match {
      case Vector() => part
      case matches =>
        val (starts, ends, delims) = matches.map { m =>
          val newDelim = m.group(1).toInt.toChar.toString
          (m.end, m.start, newDelim)
        }.unzip3
        makeNewPart(part, starts, ends, delims)
    }
  }
}
