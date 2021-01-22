package dx.translator

import java.nio.file.Path
import dx.api.{DxApi, DxProject}
import dx.core.ir._
import dx.core.languages.Language.Language
import dx.translator.wdl.WdlTranslatorFactory
import dx.translator.cwl.CwlTranslatorFactory
import dx.util.{FileSourceResolver, FileUtils, Logger}

trait Translator {

  /**
    * The name of the runtime asset the compiler must bundle
    * with generated applets.
    */
  def runtimeAssetName: String

  /**
    * The executor JAR file.
    */
  def runtimeJar: String

  def apply: Bundle

  def fileResolver: FileSourceResolver

  def translateInputs(bundle: Bundle,
                      inputs: Vector[Path],
                      defaults: Option[Path],
                      project: DxProject): (Bundle, FileSourceResolver) = {
    val inputTranslator = new InputTranslator(bundle, inputs, defaults, project, fileResolver)
    inputTranslator.writeTranslatedInputs()
    (inputTranslator.bundleWithDefaults, inputTranslator.fileResolver)
  }
}

trait TranslatorFactory {
  def create(sourceFile: Path,
             language: Option[Language],
             locked: Boolean,
             defaultRuntimeAttrs: Map[String, Value],
             reorgAttrs: ReorgSettings,
             perWorkflowAttrs: Map[String, DxWorkflowAttrs],
             defaultScatterChunkSize: Int,
             fileResolver: FileSourceResolver,
             dxApi: DxApi = DxApi.get,
             logger: Logger = Logger.get): Option[Translator]
}

object TranslatorFactory {
  private val translatorFactories: Vector[TranslatorFactory] = Vector(
      WdlTranslatorFactory(),
      CwlTranslatorFactory()
  )

  def createTranslator(source: Path,
                       language: Option[Language] = None,
                       extras: Option[Extras] = None,
                       defaultScatterChunkSize: Int,
                       locked: Boolean = false,
                       reorgEnabled: Option[Boolean] = None,
                       baseFileResolver: FileSourceResolver = FileSourceResolver.get,
                       dxApi: DxApi = DxApi.get,
                       logger: Logger = Logger.get): Translator = {
    val sourceAbsPath = FileUtils.absolutePath(source)
    val fileResolver = baseFileResolver.addToLocalSearchPath(Vector(sourceAbsPath.getParent))
    // load defaults from extras
    val defaultRuntimeAttrs = extras.map(_.defaultRuntimeAttributes).getOrElse(Map.empty)
    val reorgAttrs = (extras.flatMap(_.customReorgAttributes), reorgEnabled) match {
      case (Some(attr), None)    => attr
      case (Some(attr), Some(b)) => attr.copy(enabled = b)
      case (None, Some(b))       => DefaultReorgSettings(b)
      case (None, None)          => DefaultReorgSettings(false)
    }
    val perWorkflowAttrs = extras.map(_.perWorkflowDxAttributes).getOrElse(Map.empty)
    translatorFactories
      .collectFirst { factory =>
        factory.create(sourceAbsPath,
                       language,
                       locked,
                       defaultRuntimeAttrs,
                       reorgAttrs,
                       perWorkflowAttrs,
                       defaultScatterChunkSize,
                       fileResolver,
                       dxApi,
                       logger) match {
          case Some(translator) => translator
        }
      }
      .getOrElse(
          language match {
            case Some(lang) =>
              throw new Exception(s"Language ${lang} is not supported")
            case None =>
              throw new Exception(s"Cannot determine language/version from source file ${source}")
          }
      )
  }
}
