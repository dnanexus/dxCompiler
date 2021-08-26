package dx.core.languages.wdl

import dx.core.ir.DxName

/**
  * Simplified decoding for parameter names that conform to DNAnexus character
  * restrictions, with the possible exception of containg dots.
  */
class WdlDxName(encodedPrefix: Option[String] = None,
                decodedPrefix: Option[String] = None,
                suffix: Option[String] = None)
    extends DxName(encodedPrefix, decodedPrefix, suffix) {
  assert(!encodedPrefix.exists(_.contains(DxName.Dot)),
         s"encoded prefix ${encodedPrefix.get} contains '.'")
  assert(!decodedPrefix.exists(_.endsWith(DxName.Dot)),
         s"decoded prefix ${decodedPrefix.get} ends with a '.'")

  override protected val namespaceDelim: String = DxName.Dot

  override protected def create(encodedPrefix: Option[String] = None,
                                decodedPrefix: Option[String] = None,
                                suffix: Option[String] = None): WdlDxName = {
    new WdlDxName(encodedPrefix, decodedPrefix, suffix)
  }

  override protected def encodePrefix(prefix: String): String = {
    prefix.split("\\.").toVector match {
      case Vector(_) => prefix
      case parts     =>
        // check that none of the individual words start/end with '_',
        // which will cause problems when decoding
        if (parts.dropRight(1).exists(_.endsWith("_")) || parts.drop(1).exists(_.startsWith("_"))) {
          throw new Exception(s"cannot encode prefix ${prefix} - cannot have '_' next to '.'")
        }
        parts.mkString(DxName.DotEncoded)
    }
  }

  override def decodePrefix(prefix: String): String = {
    encoded.split(DxName.DotEncoded).toVector match {
      case Vector(_) => encoded
      case parts     =>
        // check that none of the individual words start/end with '_'
        if (parts.dropRight(1).exists(_.endsWith("_")) || parts.drop(1).exists(_.startsWith("_"))) {
          throw new Exception(
              s"cannot decode value ${encoded} - more than three consecutive underscores"
          )
        }
        parts.mkString(DxName.Dot)
    }
  }
}

object WdlDxName {
  private def parse(name: String): (String, Option[String]) = {
    DxName.suffixes
      .collectFirst {
        case suffix if name.endsWith(suffix) => (name.dropRight(suffix.length), Some(suffix))
      }
      .getOrElse((name, None))
  }

  def fromEncodedParameterName(name: String): WdlDxName = {
    val (prefix, suffix) = parse(name)
    new WdlDxName(encodedPrefix = Some(prefix), suffix = suffix)
  }

  def fromDecodedParameterName(name: String): WdlDxName = {
    val (prefix, suffix) = parse(name)
    new WdlDxName(decodedPrefix = Some(prefix), suffix = suffix)
  }

  /**
    * Creates a SimpleDxName from a name that does not contain any
    * disallowed characters.
    */
  def fromRawParameterName(name: String,
                           namespace: Option[String] = None,
                           suffix: Option[String] = None): WdlDxName = {
    namespace match {
      case Some(ns) =>
        new WdlDxName(decodedPrefix = Some(s"${ns}${DxName.Dot}${name}"), suffix = suffix)
      case None =>
        new WdlDxName(encodedPrefix = Some(name), decodedPrefix = Some(name), suffix = suffix)
    }
  }
}
