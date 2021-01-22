package dx.translator.cwl

import dx.api.DxApi
import dx.core.ir.{Bundle, Value}
import dx.core.languages.Language
import dx.core.languages.Language.Language
import dx.core.languages.cwl.{CwlUtils, DxHintSchema}
import dx.cwl.{CommandLineTool, CwlRecord, Parser, RequirementUtils}
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

  override val runtimeAssetName: String = "dxCWLrt"

  override val runtimeJar: String = "dxExecutorCwl.jar"

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
    val primaryCallable = callables.head
    val irTypeAliases = typeAliases.collect {
      case (name, record: CwlRecord) => name -> CwlUtils.toIRType(record)
    }
    Bundle(Some(primaryCallable),
           Map(primaryCallable.name -> primaryCallable),
           Vector(primaryCallable.name),
           irTypeAliases)
  }
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
    lazy val basePath = fileResolver.localSearchPath match {
      case Vector()     => sourceFile.toAbsolutePath.getParent
      case Vector(path) => path
      case v =>
        logger.warning(
            s"CWL parser can only use a single import directory; ignoring ${v.tail.mkString(",")}"
        )
        v.head
    }
    lazy val parser =
      Parser.create(baseUri = Some(basePath.toUri.toString), hintSchemas = Vector(DxHintSchema))
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
    } else if (!parser.canParse(sourceFile)) {
      // otherwise make sure the file is parseable as CWL
      return None
    }
    val tool = parser.parseFile(sourceFile) match {
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
