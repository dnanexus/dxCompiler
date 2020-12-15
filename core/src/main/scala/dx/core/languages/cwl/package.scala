package dx.core.languages.cwl

import dx.core.ir.DocumentSource
import dx.util.FileUtils

import java.nio.file.Path

case class CwlDocumentSource(source: Path) extends DocumentSource {
  override def toString: String = FileUtils.readFileContent(source)
}
