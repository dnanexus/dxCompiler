package dx.translator.wdl

import java.nio.file.Path
import dx.api.{DxApi, DxProject}
import dx.core.ir._
import dx.core.ir.Type.TSchema
import dx.core.languages.Language
import dx.core.languages.wdl.{VersionSupport, WdlBundle, WdlDxName, WdlOptions, WdlUtils}
import dx.translator.{
  DxWorkflowAttrs,
  InputTranslator,
  ReorgSettings,
  Translator,
  TranslatorFactory
}
import dx.util.{FileSourceResolver, Logger}
import spray.json.{JsArray, JsObject, JsString, JsValue}
import wdlTools.syntax.NoSuchParserException
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}

/**
  * WDL input details
  * - JSON only
  * - For workflow invocation, inputs must be prefixed with workflow name
  * - Namespaced inputs are allowed for unlocked workflows, e.g. `wf1.wf2.mytask.foo`
  * - Runtime and hints can be overridden, e.g. `"wf.mytask.runtime.cpu": 2`
  */
case class WdlInputTranslator(bundle: Bundle,
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
                            complexPathValues = false,
                            ignoreUnusedInputs = false,
                            WdlDxName,
                            baseFileResolver,
                            dxApi,
                            logger) {

  override protected def translateJsInput(jsv: JsValue, t: Type): JsValue = {
    (t, jsv) match {
      case (pairType: TSchema, JsArray(pair))
          if WdlUtils.isPairSchema(pairType) && pair.size == 2 =>
        // pair represented as [left, right]
        JsObject(WdlUtils.PairLeftKey -> pair(0), WdlUtils.PairRightKey -> pair(1))
      case (mapType: TSchema, JsObject(fields))
          if WdlUtils.isMapSchema(mapType) && !WdlUtils.isMapValue(fields) =>
        // map represented as a JSON object (i.e. has keys that are coercible from String)
        val (keys, values) = fields.map {
          case (key, value) => (JsString(key), value)
        }.unzip
        JsObject(WdlUtils.MapKeysKey -> JsArray(keys.toVector),
                 WdlUtils.MapValuesKey -> JsArray(values.toVector))
      case _ => jsv
    }
  }
}

/**
  * Compiles WDL to IR.
  * TODO: remove limitation that two callables cannot have the same name
  * TODO: rewrite sortByDependencies using a graph data structure
  */
case class WdlTranslator(doc: TAT.Document,
                         typeAliases: Map[String, WdlTypes.T_Struct],
                         locked: Boolean,
                         defaultRuntimeAttrs: Map[String, Value],
                         reorgAttrs: ReorgSettings,
                         perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                         defaultScatterChunkSize: Int,
                         useManifests: Boolean,
                         versionSupport: VersionSupport,
                         fileResolver: FileSourceResolver = FileSourceResolver.get,
                         dxApi: DxApi = DxApi.get,
                         logger: Logger = Logger.get)
    extends Translator {

  override val runtimeAssetName: String = "dxWDLrt"

  override val runtimeJar: String = "dxExecutorWdl.jar"

  override val complexPathValues: Boolean = false

  override lazy val apply: Bundle = {
    val wdlBundle: WdlBundle = WdlBundle.create(doc)
    val callableTranslator = CallableTranslator(
        wdlBundle,
        typeAliases,
        locked,
        defaultRuntimeAttrs,
        reorgAttrs,
        perWorkflowAttrs,
        defaultScatterChunkSize,
        useManifests,
        versionSupport,
        dxApi,
        fileResolver,
        logger
    )
    // sort callables by dependencies
    val logger2 = logger.withIncTraceIndent()
    val depOrder: Vector[TAT.Callable] = wdlBundle.sortByDependencies(logger2)
    if (logger2.isVerbose) {
      logger2.trace(s"all tasks: ${wdlBundle.tasks.keySet}")
      logger2.trace(s"all callables in dependency order: ${depOrder.map(_.name)}")
    }
    // translate each callable in order
    val (allCallables, sortedCallables) =
      depOrder.foldLeft((Map.empty[String, Callable], Vector.empty[Callable])) {
        case ((allCallables, sortedCallables), callable) =>
          val translatedCallables = callableTranslator.translateCallable(callable, allCallables)
          (
              allCallables ++ translatedCallables.map(c => c.name -> c).toMap,
              sortedCallables ++ translatedCallables
          )
      }
    val allCallablesSortedNames = sortedCallables.map(_.name).distinct
    val primaryCallable = wdlBundle.primaryCallable.map { callable =>
      allCallables(WdlUtils.getUnqualifiedName(callable.name))
    }
    if (logger2.isVerbose) {
      logger2.trace(s"allCallables: ${allCallables.keys}")
      logger2.trace(s"allCallablesSorted: ${allCallablesSortedNames}")
    }
    val irTypeAliases = typeAliases.map {
      case (name, struct: WdlTypes.T_Struct) => name -> WdlUtils.toIRType(struct)
    }
    Bundle(primaryCallable, allCallables, allCallablesSortedNames, irTypeAliases)
  }

  override protected def createInputTranslator(bundle: Bundle,
                                               inputs: Vector[Path],
                                               defaults: Option[Path],
                                               project: DxProject): InputTranslator = {
    WdlInputTranslator(bundle, inputs, defaults, project, useManifests, fileResolver)
  }
}

case class WdlTranslatorFactory(wdlOptions: WdlOptions = WdlOptions.default)
    extends TranslatorFactory {
  override def create(sourceFile: Path,
                      language: Option[Language.Language],
                      locked: Boolean,
                      defaultRuntimeAttrs: Map[String, Value],
                      reorgAttrs: ReorgSettings,
                      perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                      defaultScatterChunkSize: Int,
                      useManifests: Boolean,
                      fileResolver: FileSourceResolver,
                      dxApi: DxApi = DxApi.get,
                      logger: Logger = Logger.get): Option[WdlTranslator] = {
    val (doc, typeAliases, versionSupport) =
      try {
        VersionSupport.fromSourceFile(sourceFile, wdlOptions, fileResolver, dxApi, logger)
      } catch {
        // If the exception is because this is not a WDL document, return
        // None so other translators will have a chance to try to parse it,
        // otherwise it is a WDL document with a syntax or type error and
        // we let those exceptions surface.
        case _: NoSuchParserException => return None
      }
    if (!language.forall(Language.toWdlVersion(_) == doc.version.value)) {
      throw new Exception(s"WDL document ${sourceFile} is not version ${language.get}")
    }
    Some(
        WdlTranslator(
            doc,
            typeAliases.toMap,
            locked,
            defaultRuntimeAttrs,
            reorgAttrs,
            perWorkflowAttrs,
            defaultScatterChunkSize,
            useManifests,
            versionSupport,
            fileResolver,
            dxApi,
            logger
        )
    )
  }
}
