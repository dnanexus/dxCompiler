package dx.executor.wdl

import wdlTools.eval.{Hints, Meta, WdlValues}
import wdlTools.eval.WdlValues.{V, V_Object}
import wdlTools.syntax.WdlVersion
import wdlTools.types.{TypedAbstractSyntax => TAT}

case class MetaResolver(wdlVersion: WdlVersion,
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

  def getHint(name: String):
}
