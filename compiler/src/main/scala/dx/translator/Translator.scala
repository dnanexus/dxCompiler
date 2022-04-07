package dx.translator

import java.nio.file.Path
import dx.api.{DxApi, DxProject}
import dx.core.ir._
import dx.core.languages.Language.Language
import dx.core.languages.wdl.WdlOptions
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

  def fileResolver: FileSourceResolver

  def complexPathValues: Boolean

  def apply: Bundle

  protected def createInputTranslator(bundle: Bundle,
                                      inputs: Vector[Path],
                                      defaults: Option[Path],
                                      project: DxProject): InputTranslator

  def translateInputs(bundle: Bundle,
                      inputs: Vector[Path],
                      defaults: Option[Path],
                      project: DxProject): (Bundle, FileSourceResolver) = {
    val inputTranslator = createInputTranslator(bundle, inputs, defaults, project)
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
             useManifests: Boolean,
             instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
             fileResolver: FileSourceResolver,
             dxApi: DxApi = DxApi.get,
             logger: Logger = Logger.get): Option[Translator]
}

object TranslatorFactory {
  def createTranslator(source: Path,
                       language: Option[Language] = None,
                       wdlOptions: WdlOptions = WdlOptions.default,
                       extras: Option[Extras] = None,
                       defaultScatterChunkSize: Int,
                       locked: Boolean = false,
                       reorgEnabled: Option[Boolean] = None,
                       useManifests: Boolean,
                       instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
                       baseFileResolver: FileSourceResolver = FileSourceResolver.get,
                       dxApi: DxApi = DxApi.get,
                       logger: Logger = Logger.get): Translator = {
    val sourceAbsPath = FileUtils.absolutePath(source)
    val fileResolver = baseFileResolver.addToLocalSearchPath(Vector(sourceAbsPath.getParent))
    // load defaults from extras
    val defaultRuntimeAttrs = extras.flatMap(_.defaultRuntimeAttributes).getOrElse(Map.empty)
    val reorgAttrs = (extras.flatMap(_.customReorgAttributes), reorgEnabled) match {
      case (Some(attr), None)    => attr
      case (Some(attr), Some(b)) => attr.copy(enabled = b)
      case (None, Some(b))       => DefaultReorgSettings(b)
      case (None, None)          => DefaultReorgSettings(false)
    }
    val perWorkflowAttrs = extras.flatMap(_.perWorkflowDxAttributes).getOrElse(Map.empty)
    val translatorFactories = Vector(
        (wdlOptions),
        CwlTranslatorFactory()
    )
    translatorFactories
      .collectFirst { factory =>
        factory.create(
            sourceAbsPath,
            language,
            locked,
            defaultRuntimeAttrs,
            reorgAttrs,
            perWorkflowAttrs,
            defaultScatterChunkSize,
            useManifests,
            instanceTypeSelection,
            fileResolver,
            dxApi,
            logger
        ) match {
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
