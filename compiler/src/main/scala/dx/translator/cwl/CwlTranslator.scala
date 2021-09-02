package dx.translator.cwl

import dx.api.{DxApi, DxProject}
import dx.core.ir.{Bundle, InstanceTypeSelection, Type, Value}
import dx.core.languages.Language
import dx.core.languages.Language.Language
import dx.core.languages.cwl.{CwlUtils, DxHintSchema}
import dx.cwl.{CommandLineTool, CwlRecord, HintUtils, Parser}
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
                            baseFileResolver,
                            dxApi,
                            logger) {

  override protected def translateJsInput(jsv: JsValue, t: Type): JsValue = {
    (t, jsv) match {
      case (Type.TFile, JsObject(fields)) if fields.get("class").contains(JsString("File")) =>
        fields.get("location") match {
          case Some(uri: JsString) if uri.value.startsWith("dx://") => uri
          case Some(obj) =>
            try {
              obj
            } catch {
              case _: Throwable =>
                throw new Exception(s"not a valid DNAnexus link object: ${obj}")
            }
          case _ =>
            throw new Exception(
                s"dxCompiler can only translate files represented as dx:// URIs or as link objects, not ${fields}"
            )
        }
      case _ => jsv
    }
  }
}

case class CwlTranslator(tool: CommandLineTool,
                         sourceFile: Path,
                         locked: Boolean,
                         defaultRuntimeAttrs: Map[String, Value],
                         reorgAttrs: ReorgSettings,
                         perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                         defaultScatterChunkSize: Int,
                         useManifests: Boolean,
                         instanceTypeResolution: InstanceTypeSelection.InstanceTypeSelection,
                         fileResolver: FileSourceResolver = FileSourceResolver.get,
                         dxApi: DxApi = DxApi.get,
                         logger: Logger = Logger.get)
    extends Translator {

  override val runtimeAssetName: String = "dxCWLrt"

  override val runtimeJar: String = "dxExecutorCwl.jar"

  private lazy val typeAliases = HintUtils.getSchemaDefs(tool.requirements)

  override lazy val apply: Bundle = {
    val callableTranslator = ProcessTranslator(
        typeAliases,
        locked,
        defaultRuntimeAttrs,
        reorgAttrs,
        perWorkflowAttrs,
        defaultScatterChunkSize,
        instanceTypeResolution,
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
                      instanceTypeResolution: InstanceTypeSelection.InstanceTypeSelection,
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
    } else {
      // otherwise make sure the file is parseable as CWL
      parser.detectVersionAndClass(sourceFile) match {
        case Some((version, _)) if Language.parse(version) == Language.CwlV1_2 => ()
        case _ =>
          return None
      }
    }
    val tool = parser.parseFile(sourceFile) match {
      case tool: CommandLineTool => tool
      case _ =>
        throw new Exception(s"Not a command line tool ${sourceFile}")
    }
    Some(
        CwlTranslator(
            tool,
            sourceFile,
            locked,
            defaultRuntimeAttrs,
            reorgAttrs,
            perWorkflowAttrs,
            defaultScatterChunkSize,
            useManifests,
            instanceTypeResolution,
            fileResolver,
            dxApi,
            logger
        )
    )
  }
}
