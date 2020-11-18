package dx.translator.cwl

import java.nio.file.Path

import dx.cwl.{CommandLineTool, Parser, Process}
import dx.api.{DxApi, DxProject}
import dx.core.ir.{Application, Bundle, Callable, CallableAttribute, ExecutableKindApplet, Parameter, RuntimeRequirement, Type, Value}
import dx.core.languages.Language.Language
import dx.translator.{ReorgSettings, RunSpec, Translator, TranslatorFactory}
import wdlTools.util.{FileSourceResolver, Logger}
import dx.core.languages.cwl.CwlUtils
import dx.core.ir.Bundle
import dx.translator.cwl.CallableTranslator

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
    // doc.requirements.collect -> schemadefrequirements -> typealiases
    val app = callableTranslator.translateCallable(cwlBundle.primaryCallable, cwlBundle.sourceFile.getFileName.toString)
    val allCallables: Map[String, Callable] = Map((app.name, app)) // FIXME: works only for CommandLineTool
    Bundle(Option(app), allCallables, dependencies, Map.empty[String, Type])
  }

  // TODO: not needed for IR
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
    println(Parser.canParse(sourceFile))
    val doc = Parser.parse(sourceFile)
    Option(CwlTranslator(doc, sourceFile, locked, defaultRuntimeAttrs, reorgAttrs, fileResolver, dxApi))

  }


}