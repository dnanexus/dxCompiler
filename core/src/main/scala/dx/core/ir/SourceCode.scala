package dx.core.ir

/**
  * Abstraction of the source code associated with IR.
  */
trait SourceCode {

  /**
    * The workflow lanaguage of the source code.
    */
  def language: String

  /**
    * Specifies a target to execute within the source.
    * @return
    */
  def target: Option[String] = None

  /**
    * Generates the source code as a String.
    * @return
    */
  def toString: String
}
