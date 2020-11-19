/**
  * Wrappers around WDL-specific document elements, used by the Compiler
  * when generating apps/workflows.
  */
package dx.translator.cwl
import java.nio.file.Path

import dx.core.ir.DocumentSource

import scala.io.Source


case class CwlDocumentSource(sourceFile: Path)
  extends DocumentSource {

  override def toString: String = Source.fromFile(sourceFile.toString).getLines.mkString("\n") // FIXME return original YML file
}

