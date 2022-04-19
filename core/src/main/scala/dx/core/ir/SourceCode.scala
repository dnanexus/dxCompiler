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

  /**
    * Returns the source code for the Document in the doc.contents.
    * This is a hack for APPS-994 to compare Frag/Block parsed source code for applet reuse.
    * Parsed source is not in a WDL format therefore toString method does not work here
    * Suggested a permanent solution to create Frags/Blocks along with detection of Tasks and parse a .wdl correctly.
    * See the APPS-994 for comments @Gvaihir
    */
  def getDocContents: String

}
