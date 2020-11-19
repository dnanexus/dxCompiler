package dx.translator.cwl

import java.nio.file.Path

import dx.cwl.{Parser, Process}
import dx.api.{DxApi, DxProject}
import dx.core.ir.{Callable, Type, Value}
import dx.translator.{ReorgSettings, Translator, TranslatorFactory}
import dx.util.{FileSourceResolver, Logger}
import dx.core.ir.Bundle
import dx.core.languages.Language.Language

case class CwlTranslator(doc: Process,
                         sourceFile: Path,
                         locked: Boolean,
                         defaultRuntimeAttrs: Map[String, Value],
                         reorgAttrs: ReorgSettings,
                         fileResolver: FileSourceResolver = FileSourceResolver.get,
                         dxApi: DxApi = DxApi.get)
  extends Translator {

  override lazy val apply: Bundle = {

    val dependencies = Vector.empty[String] // There are no dependencies in cwl like there were in wdl
    val cwlBundle = CwlBundle.create(doc, sourceFile)
    val callableTranslator = CallableTranslator(
      cwlBundle,
      locked,
      defaultRuntimeAttrs,
      reorgAttrs,
      dxApi
    )
    val app = callableTranslator.translateCallable(cwlBundle.primaryCallable, cwlBundle.sourceFile)
    val allCallables: Map[String, Callable] = Map((app.name, app)) // TODO: workflows can have more callables
    Bundle(Option(app), allCallables, dependencies, Map.empty[String, Type])
  }

  def translateInputs(bundle: Bundle,
                      inputs: Vector[Path],
                      defaults: Option[Path],
                      project: DxProject): (Bundle, FileSourceResolver) = {
    (bundle, FileSourceResolver.get)
  }
}


case class CwlTranslatorFactory()
  extends TranslatorFactory {

  def create(sourceFile: Path,
             language: Option[Language],
             locked: Boolean,
             defaultRuntimeAttrs: Map[String, Value],
             reorgAttrs: ReorgSettings,
             fileResolver: FileSourceResolver,
             dxApi: DxApi = DxApi.get,
             logger: Logger = Logger.get): Option[Translator] = {
    val doc = Parser.parse(sourceFile)
    Option(CwlTranslator(doc, sourceFile, locked, defaultRuntimeAttrs, reorgAttrs, fileResolver, dxApi))
  }
}