package dx.translator.cwl

import dx.api.{DxApi, DxProject}
import dx.core.ir.{Bundle, Callable, Type, Value, ValueSerde}
import dx.core.languages.Language
import dx.core.languages.Language.Language
import dx.core.languages.cwl.{CwlBundle, CwlUtils, DxHintSchema}
import dx.cwl.{
  CommandLineTool,
  CwlDirectory,
  CwlFile,
  CwlRecord,
  DirectoryValue,
  FileValue,
  Parser,
  Process,
  Workflow
}
import dx.translator.{
  DxWorkflowAttrs,
  InputTranslator,
  ReorgSettings,
  Translator,
  TranslatorFactory
}
import dx.util.{FileSourceResolver, Logger}
import org.w3id.cwl.cwl1_2.CWLVersion
import spray.json._

import java.nio.file.Path

case class CwlInputTranslator(bundle: Bundle,
                              inputs: Vector[Path],
                              defaults: Option[Path],
                              project: DxProject,
                              useManifests: Boolean,
                              baseFileResolver: FileSourceResolver = FileSourceResolver.get,
                              dxApi: DxApi = DxApi.get,
                              logger: Logger = Logger.get)
    extends InputTranslator(bundle,
                            inputs,
                            defaults,
                            project,
                            useManifests,
                            complexPathValues = true,
                            ignoreUnusedInputs = true,
                            baseFileResolver,
                            dxApi,
                            logger) {

  override protected def translateJsInput(jsv: JsValue, t: Type): JsValue = {
    (t, jsv) match {
      case (Type.TFile, obj: JsObject) if obj.fields.get("class").contains(JsString("File")) =>
        val (_, cwlValue) = CwlUtils.toIRValue(FileValue.deserialize(obj), CwlFile)
        ValueSerde.serialize(cwlValue, fileResolver = Some(baseFileResolver), pathsAsObjects = true)
      case (Type.TDirectory, obj: JsObject)
          if obj.fields.get("class").contains(JsString("Directory")) =>
        val (_, cwlValue) = CwlUtils.toIRValue(DirectoryValue.deserialize(obj), CwlDirectory)
        ValueSerde.serialize(cwlValue, fileResolver = Some(baseFileResolver), pathsAsObjects = true)
      case _ => jsv
    }
  }
}

case class CwlTranslator(process: Process,
                         sourceFile: Path,
                         locked: Boolean,
                         defaultRuntimeAttrs: Map[String, Value],
                         reorgAttrs: ReorgSettings,
                         perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                         defaultScatterChunkSize: Int,
                         useManifests: Boolean,
                         fileResolver: FileSourceResolver = FileSourceResolver.get,
                         dxApi: DxApi = DxApi.get,
                         logger: Logger = Logger.get)
    extends Translator {

  override val runtimeAssetName: String = "dxCWLrt"

  override val runtimeJar: String = "dxExecutorCwl.jar"

  override val complexPathValues: Boolean = true

  override lazy val apply: Bundle = {
    val cwlBundle: CwlBundle = CwlBundle.create(process)
    val callableTranslator = ProcessTranslator(
        cwlBundle,
        locked,
        defaultRuntimeAttrs,
        reorgAttrs,
        perWorkflowAttrs,
        defaultScatterChunkSize,
        useManifests,
        dxApi,
        fileResolver,
        logger
    )
    // sort callables by dependencies
    val logger2 = logger.withIncTraceIndent()
    val depOrder: Vector[Process] = cwlBundle.sortByDependencies
    if (logger2.isVerbose) {
      logger2.trace(s"all tasks: ${cwlBundle.tools.keySet}")
      logger2.trace(s"all processes in dependency order: ${depOrder.map(_.name)}")
    }
    // translate processes
    val (allCallables, sortedCallables) =
      depOrder.foldLeft((Map.empty[String, Callable], Vector.empty[Callable])) {
        case ((allCallables, sortedCallables), callable) =>
          val isPrimary = callable.name == cwlBundle.primaryProcess.name
          val translatedCallables =
            callableTranslator.translateProcess(callable, allCallables, isPrimary = isPrimary)
          (
              allCallables ++ translatedCallables.map(c => c.name -> c).toMap,
              sortedCallables ++ translatedCallables
          )
      }
    val allCallablesSortedNames = sortedCallables.map(_.name).distinct
    val primaryCallable = allCallables(cwlBundle.primaryProcess.name)
    if (logger2.isVerbose) {
      logger2.trace(s"allCallables: ${allCallables.keys}")
      logger2.trace(s"allCallablesSorted: ${allCallablesSortedNames}")
    }
    val irTypeAliases = cwlBundle.typeAliases.collect {
      case (name, record: CwlRecord) => name -> CwlUtils.toIRType(record)
    }
    Bundle(Some(primaryCallable), allCallables, allCallablesSortedNames, irTypeAliases)
  }

  override protected def createInputTranslator(bundle: Bundle,
                                               inputs: Vector[Path],
                                               defaults: Option[Path],
                                               project: DxProject): InputTranslator = {
    CwlInputTranslator(bundle, inputs, defaults, project, useManifests, fileResolver)
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
                      useManifests: Boolean,
                      fileResolver: FileSourceResolver,
                      dxApi: DxApi,
                      logger: Logger): Option[Translator] = {
    // TODO: we need to require that the source file be "packed" before compiling, because
    //  we cannot include auxiliary files (e.g. a JavaScript or YAML import) with the CWL.
    //  Then we shouldn't use a base URI and instead let parsing errors due to unsatisfied
    //  imports (which shouldn't happen) bubble up. We should also print a warning if the
    //  user tries to specify any import directories for CWL.
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
      Parser.create(baseUri = Some(basePath.toUri), hintSchemas = Vector(DxHintSchema))
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
    } else {
      // otherwise make sure the file is parseable as CWL
      parser.detectVersionAndClass(sourceFile) match {
        case Some((version, _)) if Language.parse(version) == Language.CwlV1_2 => ()
        case _ =>
          return None
      }
    }
    val process = parser.parseFile(sourceFile) match {
      case tool: CommandLineTool => tool
      case wf: Workflow          => wf
      case _ =>
        throw new Exception(s"Not a command line tool or workflow: ${sourceFile}")
    }
    Some(
        CwlTranslator(
            process,
            sourceFile,
            locked,
            defaultRuntimeAttrs,
            reorgAttrs,
            perWorkflowAttrs,
            defaultScatterChunkSize,
            useManifests,
            fileResolver,
            dxApi,
            logger
        )
    )
  }
}
