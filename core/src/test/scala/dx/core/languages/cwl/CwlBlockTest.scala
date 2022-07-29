package dx.core.languages.cwl

import dx.core.io.DxWorkerPaths
import dx.core.ir.Workflow
import dx.core.languages.Language
import dx.cwl.{Parser, ParserResult, Process, Workflow}
import dx.util.PosixPath
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.{Files, Path, Paths}

class CwlBlockTest extends AnyFlatSpec with Matchers {

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  private def getWfBlocks(processName: String): Vector[CwlBlock] = {
    val cwlFile: Path = pathFromBasename("cwl", s"${processName}.cwl.json")
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()
    val workerPaths = DxWorkerPaths(PosixPath(jobRootDir.toString))
    workerPaths.createCleanDirs()

    val parser = Parser.create(hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClassFromFile(cwlFile) match {
      case (version, _) if Language.parse(version) == Language.CwlV1_2 => ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL document of a supported version
               |${cwlFile}""".stripMargin
        )
    }
    val process = parser.parseFile(cwlFile) match {
      case ParserResult(Some(process: Process), _, _, _) => process
      case other =>
        throw new Exception(s"expected CWL document to contain a CommandLineTool, not ${other}")
    }
    process match {
      case wf: Workflow => CwlBlock.createBlocks(wf)
      case _            => throw new Exception(s"Unexpected type of process ${process}")
    }
  }

}
