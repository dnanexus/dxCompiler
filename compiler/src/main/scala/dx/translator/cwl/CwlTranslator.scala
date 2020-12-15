package dx.translator.cwl

import dx.api.{DxApi, DxProject}
import dx.core.ir.{Bundle, Value}
import dx.core.languages.Language
import dx.core.languages.Language.Language
import dx.core.languages.cwl.CwlUtils
import dx.cwl.{CommandLineTool, CwlRecord, Parser, Requirement, RequirementUtils}
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

  private lazy val typeAliases = RequirementUtils.getSchemaDefs(tool.requirements)

  override lazy val apply: Bundle = {
    val callableTranslator = ProcessTranslator(
        typeAliases,
        locked,
        defaultRuntimeAttrs,
        reorgAttrs,
        perWorkflowAttrs,
        defaultScatterChunkSize,
        dxApi,
        fileResolver,
        logger
    )
    val callables = callableTranslator.translateProcess(tool)
    assert(callables.size == 1)
    val irTypeAliases = typeAliases.collect {
      case (name, record: CwlRecord) => name -> CwlUtils.toIRType(record)
    }
    Bundle(Some(callables.head), Map.empty, Vector.empty, irTypeAliases)
  }

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
    val tool = Parser.parseFile(sourceFile, hintSchemas = DxHintSchema) match {
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
