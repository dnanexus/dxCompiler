package dx.translator.cwl

import dx.api.{DxApi, DxProject}
import dx.core.ir.{Bundle, Value}
import dx.core.languages.Language
import dx.core.languages.Language.Language
import dx.cwl.{CommandLineTool, Parser}
import dx.translator.{DxWorkflowAttrs, ReorgSettings, Translator, TranslatorFactory}
import dx.util.{FileSourceResolver, Logger}
import org.w3id.cwl.cwl1_2.CWLVersion

import java.nio.file.Path

case class CwlTranslator(tool: CommandLineTool,
                         sourceFile: Path,
                         locked: Boolean,
                         defaultRuntimeAttrs: Map[String, Value],
                         reorgAttrs: ReorgSettings,
                         perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                         defaultScatterChunkSize: Int,
                         fileResolver: FileSourceResolver = FileSourceResolver.get,
                         dxApi: DxApi = DxApi.get,
                         logger: Logger = Logger.get)
    extends Translator {

  /**
    * The name of the runtime asset the compiler must bundle
    * with generated applets.
    */
  override def runtimeAssetName: String = "dxCWLrt"

  override lazy val apply: Bundle = {}

  override def translateInputs(bundle: Bundle,
                               inputs: Vector[Path],
                               defaults: Option[Path],
                               project: DxProject): (Bundle, FileSourceResolver) = {}
}

case class CwlTranslatorFactory() extends TranslatorFactory {
  override def create(sourceFile: Path,
                      language: Option[Language],
                      locked: Boolean,
                      defaultRuntimeAttrs: Map[String, Value],
                      reorgAttrs: ReorgSettings,
                      perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                      defaultScatterChunkSize: Int,
                      fileResolver: FileSourceResolver,
                      dxApi: DxApi,
                      logger: Logger): Option[Translator] = {
    if (language.isDefined) {
      // if language is specified, make sure it is CWL 1.2
      val ver =
        try {
          Language.toCwlVersion(language.get)
        } catch {
          case _: Throwable =>
            return None
        }
      if (ver != CWLVersion.V1_2) {
        throw new Exception(s"dxCompiler does not support CWL version ${ver}")
      }
    } else if (!Parser.canParse(sourceFile)) {
      // otherwise make sure the file is parseable as CWL
      return None
    }
    val tool = Parser.parse(sourceFile) match {
      case tool: CommandLineTool => tool
      case _ =>
        throw new Exception(s"Not a command line tool ${sourceFile}")
    }
    Some(
        CwlTranslator(tool,
                      sourceFile,
                      locked,
                      defaultRuntimeAttrs,
                      reorgAttrs,
                      perWorkflowAttrs,
                      defaultScatterChunkSize,
                      fileResolver,
                      dxApi,
                      logger)
    )
  }
}
