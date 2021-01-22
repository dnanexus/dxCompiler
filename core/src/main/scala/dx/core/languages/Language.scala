package dx.core.languages

import wdlTools.syntax.WdlVersion
import dx.util.Enum
import org.w3id.cwl.cwl1_2.CWLVersion

object Language extends Enum {
  type Language = Value
  val WdlVDraft2, WdlV1_0, WdlV2_0, CwlV1_2 = Value
  val WdlDefault: Language = WdlV1_0
  val CwlDefault: Language = CwlV1_2

  def toWdlVersion(value: Language): WdlVersion = {
    value match {
      case WdlVDraft2 => WdlVersion.Draft_2
      case WdlV1_0    => WdlVersion.V1
      case WdlV2_0    => WdlVersion.V2
      case other      => throw new Exception(s"${other} is not a WDL version")
    }
  }

  def fromWdlVersion(version: WdlVersion): Value = {
    version match {
      case WdlVersion.Draft_2 => Language.WdlVDraft2
      case WdlVersion.V1      => Language.WdlV1_0
      case WdlVersion.V2      => Language.WdlV2_0
      case other              => throw new Exception(s"Unsupported WDL version ${other}")
    }
  }

  def toCwlVersion(value: Language): CWLVersion = {
    value match {
      case CwlV1_2 => CWLVersion.V1_2
      case other   => throw new Exception(s"${other} is not a CWL version")
    }
  }

  def fromCwlVersion(version: CWLVersion): Value = {
    version match {
      case CWLVersion.V1_2 => Language.CwlV1_2
      case other           => throw new Exception(s"Unsupported CWL version ${other}")
    }
  }

  private def normalizeVersion(version: String): String = {
    version.trim.toLowerCase.replaceAll("[._-]", "")
  }

  private val languageRegexp = s"(cwl|wdl)[_ ]?(.*)".r

  /**
    * Parses a version string.
    * @param version The version specifier. Should be "<language> <version>", but the language may be omitted if
    *                `hint` is specified, or if the version is unambiguous.
    * @param languageHint Hint about which language is being used, when `version`
    *                     does not include the language specifier.
    * @return
    */
  def parse(version: String, languageHint: Option[String] = None): Value = {
    val (lang: Option[String], ver: Option[String]) = version.toLowerCase match {
      case languageRegexp(l, _) if languageHint.nonEmpty && l != languageHint.get =>
        throw new Exception(s"language hints don't match ${l} ${languageHint.get}")
      case languageRegexp(l, v) if v.nonEmpty =>
        (Some(l), Some(normalizeVersion(v)))
      case languageRegexp(l, _) =>
        (Some(l), None)
      case _ if version.nonEmpty =>
        (languageHint.map(normalizeVersion), Some(normalizeVersion(version)))
      case _ =>
        (languageHint.map(normalizeVersion), None)
    }
    (lang, ver) match {
      case (Some("wdl"), None)                                      => WdlDefault
      case (Some("cwl"), None)                                      => CwlDefault
      case (None | Some("wdl"), Some("draft2"))                     => Language.WdlVDraft2
      case (None | Some("wdl"), Some("draft3" | "10" | "v10"))      => Language.WdlV1_0
      case (None | Some("wdl"), Some("development" | "20" | "v20")) => Language.WdlV2_0
      case (None | Some("cwl"), Some("12" | "v12"))                 => Language.CwlV1_2
      case other =>
        throw new Exception(s"Unrecognized/unsupported language ${other}")
    }
  }
}
