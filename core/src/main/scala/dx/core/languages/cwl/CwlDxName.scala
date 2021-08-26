package dx.core.languages.cwl

import dx.core.ir.DxName

object CwlDxName {
  val Delim = "/"

  // character sequences that may not appear in a non-encoded name
  private val illegalSeqsRegexp = "__|\\s+".r
  // characters that need to be encoded
  private val encodeCharsRegexp = "([^a-zA-Z0-9_])".r
  // character sequences that need to be decoded
  private val decodeSeqsRegexp = "(___)|__(\\d+)__".r

  /**
    * Creates a ComplexDxName from an encoded parameter name that may
    * have a suffix.
    */
  def fromEncodedParameterName(name: String): CwlDxName = {
    val (prefix, suffix) = DxName.split(name)
    new CwlDxName(encodedPrefix = Some(prefix), suffix = suffix)
  }

  /**
    * Creates a ComplexDxName from a decoded parameter name that may
    * have a suffix.
    */
  def fromDecodedParameterName(name: String): CwlDxName = {
    val (prefix, suffix) = DxName.split(name)
    new CwlDxName(decodedPrefix = Some(prefix), suffix = suffix)
  }

  /**
    * Creates a ComplexDxName from a parameter name coming directly from a
    * language AST, with optional namespace.
    */
  def fromRawParameterName(name: String,
                           namespace: Option[String] = None,
                           suffix: Option[String] = None): CwlDxName = {
    new CwlDxName(
        decodedPrefix = Some(namespace.map(prefix => s"${prefix}${Delim}${name}").getOrElse(name)),
        suffix = suffix
    )
  }
}

class CwlDxName(encodedPrefix: Option[String] = None,
                decodedPrefix: Option[String] = None,
                suffix: Option[String] = None)
    extends DxName(encodedPrefix, decodedPrefix, suffix) {
  assert(!encodedPrefix.exists(_.contains(DxName.Dot)),
         s"encoded prefix ${encodedPrefix.get} contains '.'")
  assert(!decodedPrefix.exists(_.endsWith(DxName.Dot)),
         s"decoded prefix ${decodedPrefix.get} ends with a '.'")
  encodedPrefix.flatMap(CwlDxName.encodeCharsRegexp.findFirstIn).map { c =>
    throw new Exception(s"encoded prefix ${encodedPrefix.get} contains illegal character '${c}'")
  }
  decodedPrefix.flatMap(CwlDxName.illegalSeqsRegexp.findFirstIn).map { c =>
    throw new Exception(s"parameter name ${decoded} contains illegal character sequence '${c}'")
  }

  override protected def namespaceDelim: String = CwlDxName.Delim

  override protected def create(encodedPrefix: Option[String] = None,
                                decodedPrefix: Option[String] = None,
                                suffix: Option[String] = None): CwlDxName = {
    new CwlDxName(encodedPrefix, decodedPrefix, suffix)
  }

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
    * Encodes disallowed characters in a parameter name. DNAnexus only allows
    * [a-zA-Z0-9_] in parameter names. We also disallow:
    *   * whitespace - "foo bar"
    *   * consecutive underscores - "foo\_\_bar"
    *   * name segments that begin or end in underscore - "foo_._bar"
    * @param prefix the decoded prefix
    * @return the encoded parameter name, with "." replaced with "\_\_\_"
    *         and any other illegal character replaced with "\_\_ord(x)\_\_",
    *         where `ord(x)` is the ordinal value of the illegal character.
    */
  override protected def encodePrefix(prefix: String): String = {
    CwlDxName.encodeCharsRegexp.findAllMatchIn(prefix).toVector match {
      case Vector() => prefix
      case matches =>
        val (starts, ends, delims) = matches.map { m =>
          val newDelim = m.group(1) match {
            case DxName.Dot => DxName.DotEncoded
            case c if c.length() == 1 =>
              s"__${c.charAt(0).toInt}__"
            case other =>
              throw new Exception(s"unexpected multi-character delimiter ${other}")
          }
          (m.end, m.start, newDelim)
        }.unzip3
        makeNewName(prefix, starts, ends, delims)
    }
  }

  override protected def decodePrefix(prefix: String): String = {
    CwlDxName.decodeSeqsRegexp.findAllMatchIn(prefix).toVector match {
      case Vector() => prefix
      case matches =>
        val (starts, ends, delims) = matches.map { m =>
          val newDelim = (m.group(1), m.group(2)) match {
            case (DxName.DotEncoded, null) => DxName.Dot
            case (null, ord)               => ord.toInt.toChar.toString
            case other =>
              throw new Exception(s"unexpected delimiter ${other}")
          }
          (m.end, m.start, newDelim)
        }.unzip3
        makeNewName(prefix, starts, ends, delims)
    }
  }
}
