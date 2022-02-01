package dx.translator.cwl

import dx.api.{DxApi, DxProject}
import dx.core.ir.{Bundle, Callable, InstanceTypeSelection, Type, Value, ValueSerde}
import dx.core.languages.Language
import dx.core.languages.Language.Language
import dx.core.languages.cwl.{CwlBundle, CwlDxName, CwlUtils, DxHintSchema}
import dx.cwl.{
  CommandLineTool,
  CwlDirectory,
  CwlFile,
  CwlRecord,
  DirectoryValue,
  ExpressionTool,
  FileValue,
  Parser,
  ParserResult,
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

/**
  * CWL input details:
  * - YAML or JSON
  * - Only top-level inputs can be specified
  * - `cwl` namespace is supported - ignore any keys other than `cwl:requirements` and `cwl:hints`
  */
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
                            CwlDxName,
                            baseFileResolver,
                            dxApi,
                            logger) {

  private val CwlRegex = "^(?:.*\\.)?cwl:(.+)$".r

  override protected def splitInputs(
      rawInputs: Map[String, JsValue]
  ): (Map[String, JsValue], Map[String, JsObject]) = {
    val (values, overrides) =
      rawInputs.foldLeft(Map.empty[String, JsValue], Map.empty[String, JsValue]) {
        case ((values, overrides), (CwlRegex(key), value)) =>
          (values, overrides + (key -> value))
        case ((values, overrides), (key, value)) =>
          (values + (key -> value), overrides)
      }
    (values, if (overrides.nonEmpty) {
      val overridesObj = JsObject(overrides)
      bundle.primaryCallable match {
        case Some(c: Callable) => Map(c.name -> overridesObj)
        case _                 => bundle.allCallables.keys.map(_ -> overridesObj).toMap
      }
    } else {
      Map.empty[String, JsObject]
    })
  }

  override protected def convertRawInput(rawInput: JsValue, t: Type): JsValue = {
    (t, rawInput) match {
      case (Type.TFile, obj: JsObject) if obj.fields.get("class").contains(JsString("File")) =>
        val (_, cwlValue) = CwlUtils.toIRValue(FileValue.deserialize(obj), CwlFile)
        ValueSerde.serialize(cwlValue, fileResolver = Some(baseFileResolver), pathsAsObjects = true)
      case (m: Type.TMulti, obj: JsObject)
          if m.contains(Type.TFile) && obj.fields.get("class").contains(JsString("File")) =>
        val (_, cwlValue) = CwlUtils.toIRValue(FileValue.deserialize(obj), CwlFile)
        ValueSerde.serialize(cwlValue, fileResolver = Some(baseFileResolver), pathsAsObjects = true)
      case (Type.TDirectory, obj: JsObject)
          if obj.fields.get("class").contains(JsString("Directory")) =>
        val (_, cwlValue) = CwlUtils.toIRValue(DirectoryValue.deserialize(obj), CwlDirectory)
        ValueSerde.serialize(cwlValue, fileResolver = Some(baseFileResolver), pathsAsObjects = true)
      case (m: Type.TMulti, obj: JsObject)
          if m.contains(Type.TDirectory) &&
            obj.fields.get("class").contains(JsString("Directory")) =>
        val (_, cwlValue) = CwlUtils.toIRValue(DirectoryValue.deserialize(obj), CwlDirectory)
        ValueSerde.serialize(cwlValue, fileResolver = Some(baseFileResolver), pathsAsObjects = true)
      case _ => rawInput
    }
  }

  override protected val lockedWorkflowPrefixOptional: Boolean = true
}

case class CwlTranslator(process: Process,
                         sourceFile: Path,
                         cwlSchemas: Option[JsValue],
                         locked: Boolean,
                         defaultRuntimeAttrs: Map[String, Value],
                         reorgAttrs: ReorgSettings,
                         perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                         defaultScatterChunkSize: Int,
                         useManifests: Boolean,
                         instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
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
        cwlSchemas,
        locked,
        defaultRuntimeAttrs,
        reorgAttrs,
        perWorkflowAttrs,
        defaultScatterChunkSize,
        useManifests,
        instanceTypeSelection,
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
                      instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
                      fileResolver: FileSourceResolver,
                      dxApi: DxApi,
                      logger: Logger): Option[Translator] = {
    // TODO: we require that the source file be "packed" before compiling, because we cannot include
    //  auxiliary files (e.g. a JavaScript or YAML import) with the CWL. Then we shouldn't use a
    //  base URI and instead let parsing errors due to unsatisfied imports (which shouldn't happen)
    //  bubble up. We should also print a warning if the user tries to specify any import directories.
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
      parser.detectVersionAndClassFromFile(sourceFile) match {
        case (version, _) if Language.parse(version, Some("cwl")) == Language.CwlV1_2 => ()
        case _ =>
          return None
      }
    }
    // TODO: in the case of a document with a $graph element, if there are multiple top-level
    //  processes we first look for a process called 'main', and then we look for a single workflow.
    //  If we cannot determine the top-level process we fail with an error. We should also add an
    //  option for the user to specify the process name at compile time.
    val (process, schemas) = parser.parseFile(sourceFile) match {
      case ParserResult(Some(tool: CommandLineTool), _, _, schemas) => (tool, schemas)
      case ParserResult(Some(tool: ExpressionTool), _, _, schemas)  => (tool, schemas)
      case ParserResult(Some(wf: Workflow), _, _, schemas)          => (wf, schemas)
      case _ =>
        throw new Exception(s"Not a tool or workflow: ${sourceFile}")
    }
    Some(
        CwlTranslator(
            process,
            sourceFile,
            schemas,
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
        )
    )
  }
}
