package dx.core.ir

import spray.json.{JsNull, JsValue}

/**
  * Abstraction of the source code associated with IR.
  */
trait SourceCode {

  /**
    * The workflow lanaguage of the source code.
    */
  def language: String

  /**
    * Generates the source code as a String.
    * @return
    */
  def toString: String

  /**
    * Serializes any runtime options that are required to
    * parse the source code at execution time.
    */
  def optionsToJson: JsValue = JsNull
}
