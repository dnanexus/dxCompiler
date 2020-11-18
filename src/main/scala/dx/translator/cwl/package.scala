/**
  * Wrappers around WDL-specific document elements, used by the Compiler
  * when generating apps/workflows.
  */
package dx.translator.cwl
import dx.cwl.{Process}
import dx.core.ir.{DocumentSource}

case class CwlDocumentSource(doc: Process)
  extends DocumentSource {
  override def toString: String = "Temporary string, needs fixing" // FIXME needs to generate document.
}

