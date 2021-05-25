package dx.executor.wdl

import dx.core.languages.wdl.DxMetaHints
import wdlTools.eval.{Hints, Meta, WdlValues}
import wdlTools.eval.WdlValues.{V, V_Boolean, V_Object, V_String}
import wdlTools.syntax.WdlVersion
import wdlTools.types.WdlTypes.T_Boolean
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}

case class HintResolver(wdlVersion: WdlVersion,
                        paramterMetaSection: Option[TAT.MetaSection],
                        hintsSection: Option[TAT.MetaSection]) {
  private lazy val meta: Meta = Meta.create(wdlVersion, paramterMetaSection)
  private lazy val hints: Meta = Meta.create(wdlVersion, hintsSection)
  private lazy val hintInputs: Option[Map[String, V]] =
    hints.get(Hints.InputsKey) match {
      case Some(V_Object(inputs)) => Some(inputs)
      case _                      => None
    }
//  private lazy val hintOutputs: Option[Map[String, V]] =
//    hints.get(Hints.OutputsKey) match {
//      case Some(V_Object(outputs)) => Some(outputs)
//      case _                       => None
//    }

  def getInput(name: String): Option[WdlValues.V] = {
    meta.get(name).orElse(hintInputs.flatMap(_.get(name)))
  }

  def get(name: String, wdlTypes: Vector[WdlTypes.T]): Option[WdlValues.V] = {
    hints.get(name, wdlTypes)
  }

  lazy val globalLocalizationOptional: Boolean = {
    get(Hints.LocalizationOptionalKey, Vector(T_Boolean)) match {
      case Some(V_Boolean(localizationOptional)) => localizationOptional
      case _                                     => false
    }
  }

  def isLocalizationOptional(parameterName: String): Boolean = {
    getInput(parameterName) match {
      case Some(V_String(DxMetaHints.ParameterMetaStream)) =>
        true
      case Some(V_Object(fields)) =>
        // This enables the stream annotation in the object form of metadata value, e.g.
        // bam_file : {
        //   stream : true
        // }
        // We also support two aliases, dx_stream and localizationOptional
        fields.view
          .filterKeys(
              Set(DxMetaHints.ParameterMetaStream,
                  DxMetaHints.ParameterHintStream,
                  Hints.LocalizationOptionalKey)
          )
          .values
          .exists {
            case V_Boolean(b) => b
            case _            => false
          }
      case _ => globalLocalizationOptional
    }
  }
}
