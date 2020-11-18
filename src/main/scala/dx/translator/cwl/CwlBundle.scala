package dx.translator.cwl

import java.nio.file.Path

import dx.cwl.{CommandLineTool, Process, Workflow}
import org.w3id.cwl.cwl1_2.CWLVersion



case class CwlBundle(version: CWLVersion, primaryCallable: Process, sourceFile:Path)


object CwlBundle {

  private def bundleInfoFromDoc(doc: Process, sourceFile: Path): CwlBundle = {
    doc match {
      case t: CommandLineTool => {
        val cwlVersion = t.cwlVersion match {
          case Some(v) => v
          case None => throw new Exception("CWL version needs to be defined in source code.")
        }
        CwlBundle(cwlVersion, doc, sourceFile) }
      case _: Workflow => throw new NotImplementedError("Workflows are not supported.")
      case other => throw new NotImplementedError(s"Type ${other.getClass} is not supported/")
    }
  }


  // recurse into the imported packages
  // Check the uniqueness of tasks, Workflows, and Types
  // merge everything into one bundle.
  def create(doc: Process, sourceFile: Path): CwlBundle = {
  val cwlBundle = bundleInfoFromDoc(doc, sourceFile)
    cwlBundle
  }
}

