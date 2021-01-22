package dx.core.languages.cwl

import dx.core.ir.{CallableAttribute, ParameterAttribute}
import dx.cwl.{CwlSchema, Hint, HintSchema}

case class DxHints() extends Hint {
  def getParameterAttributes: Map[String, Vector[ParameterAttribute]] = {
    Map.empty
  }

  def getCallableAttributes: Map[String, Vector[CallableAttribute]] = {
    Map.empty
  }
}

object DxHintSchema extends HintSchema {
  override def apply(hint: Map[String, Any], schemaDefs: Map[String, CwlSchema]): Hint = {
    DxHints()
  }
}
