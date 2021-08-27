package dx.translator

import java.nio.file.{Path, Paths}
import dx.api.{DxApi, DxFile, DxFileDescCache, DxPath, DxProject}
import dx.core.Constants
import dx.core.ir.Type._
import dx.core.ir.{
  Application,
  Bundle,
  Callable,
  DxNameFactory,
  Manifest,
  ParameterLinkDeserializer,
  ParameterLinkSerializer,
  StaticInput,
  Type,
  Value,
  ValueSerde,
  Workflow
}
import dx.util.{FileSourceResolver, FileUtils, JsUtils, Logger}
import dx.util.protocols.DxFileAccessProtocol
import spray.json.{JsArray, JsNull, JsObject, JsString, JsValue}

/**
  * Tracks which keys are accessed in a map and ensures all keys are accessed exactly once.
  */
case class ExactlyOnce(name: String, fields: Map[String, JsValue], logger: Logger) {
  private var retrievedKeys: Set[String] = Set.empty

  def get(fqn: String): Option[JsValue] = {
    fields.get(fqn) match {
      case None =>
        logger.trace(s"getExactlyOnce ${fqn} => None")
        None
      case Some(v: JsValue) if retrievedKeys.contains(fqn) =>
        logger.trace(
            s"getExactlyOnce ${fqn} => Some(${v}); value already retrieved so returning None"
        )
        None
      case Some(v: JsValue) =>
        logger.trace(s"getExactlyOnce ${fqn} => Some(${v})")
        retrievedKeys += fqn
        Some(v)
    }
  }

  def checkAllUsed(): Unit = {
    val unused = fields.keySet -- retrievedKeys
    if (unused.nonEmpty) {
      throw new Exception(s"Could not map all ${name} fields. These were left: ${unused}")
    }
  }
}

abstract class InputTranslator(bundle: Bundle,
                               inputs: Vector[Path],
                               defaults: Option[Path],
                               project: DxProject,
                               useManifests: Boolean,
                               complexPathValues: Boolean,
                               ignoreUnusedInputs: Boolean,
                               dxNameFactory: DxNameFactory,
                               baseFileResolver: FileSourceResolver = FileSourceResolver.get,
                               dxApi: DxApi = DxApi.get,
                               logger: Logger = Logger.get) {

  private def loadJsonFileWithComments(path: Path): Map[String, JsValue] = {
    // skip comment lines, which start with ##
    JsUtils
      .getFields(JsUtils.jsFromFile(path))
      .view
      .filterKeys(!_.startsWith("##"))
      .toMap
  }

  private lazy val rawInputsJs =
    inputs.map(path => path -> loadJsonFileWithComments(path)).toMap
  private lazy val defaultsJs =
    defaults.map(loadJsonFileWithComments).getOrElse(Map.empty)
  private val ManifestSuffix = s".${Constants.InputManifest.decoded}"

  // All of the inputs in JSON format. Parameter names are of the form
  // <prefix>.<decoded_name>, where prefix is the task or workflow name,
  // and decoded_name is the parameter name as it appears in the source.
  private lazy val inputsJs: Map[Path, Map[String, JsValue]] = {
    rawInputsJs.map {
      case (path, jsValues) if jsValues.size == 1 && jsValues.keys.head.endsWith(ManifestSuffix) =>
        val (key, fields) = jsValues.head
        val prefix = key.dropRight(ManifestSuffix.length)
        val manifest = Manifest.parse(fields, dxNameFactory)
        path -> manifest.jsValues.map {
          case (name, value) => s"${prefix}.${name}" -> value
        }
      case (path, jsValues) => path -> jsValues
    }
  }

  /**
    * Converts a language-specific JSON value to one that can be deserialized to an IR Value.
    * This is the mechanism to support compile-time inputs/defaults that do not map implicitly
    * to an IR type.
    */
  protected def translateJsInput(jsv: JsValue, t: Type): JsValue = jsv

  /**
    * Extract Dx files from a JSON input.
    * Can be over-ridden to extract files in a language-specific way.
    * @param t the parameter type
    * @param jsv the JSON value
    * @return a Vector of JSON values that describe Dx files
    */
  private def extractDxFiles(jsv: JsValue, t: Type): Vector[JsValue] = {
    def extractFromArray(a: JsValue): Vector[JsValue] = {
      a match {
        case JsArray(paths) =>
          paths.flatMap {
            case obj: JsObject =>
              obj.fields.get("type") match {
                case Some(JsString("File"))    => extractDxFiles(obj, TFile)
                case Some(JsString("Listing")) => extractDxFiles(obj, TDirectory)
                case _                         => Vector.empty
              }
            case other => throw new Exception(s"invalid path value ${other}")
          }
        case other =>
          throw new Exception(s"invalid path array ${other}")
      }
    }

    val updatedValue = translateJsInput(jsv, t)
    (t, updatedValue) match {
      case (_, JsNull) if Type.isOptional(t) =>
        Vector.empty
      case (TOptional(inner), _) =>
        extractDxFiles(updatedValue, inner)
      case (TFile, uri: JsString) if uri.value.startsWith(DxPath.DxUriPrefix) =>
        // a dx://file-xxx or dx://project-xxx:/path/to/file URI
        Vector(uri)
      case (TFile, obj: JsObject) if DxFile.isLinkJson(obj) =>
        Vector(obj)
      case (TFile, JsObject(fields)) if !fields.contains("contents") =>
        Vector(fields("uri")) ++ fields
          .get("secondaryFiles")
          .map(extractFromArray)
          .getOrElse(Vector.empty)
      case (TDirectory, JsObject(fields)) if fields.contains("listing") =>
        fields
          .get("listing")
          .map(extractFromArray)
          .getOrElse(Vector.empty)
      case _ if Type.isPrimitive(t) => Vector.empty
      case (TArray(elementType, _), JsArray(array)) =>
        array.flatMap(element => extractDxFiles(element, elementType))
      case (TSchema(name, members), JsObject(fields)) =>
        members.flatMap {
          case (memberName, memberType) =>
            fields.get(memberName) match {
              case Some(memberValue) =>
                extractDxFiles(memberValue, memberType)
              case None if Type.isOptional(memberType) =>
                Vector.empty
              case _ =>
                throw new Exception(s"missing value for struct ${name} member ${memberName}")
            }
        }.toVector
      case (THash, JsObject(_)) =>
        // anonymous objects will never result in file-typed members, so just skip these
        Vector.empty
      case _ => Vector.empty
    }
  }

  // scan all inputs and defaults and extract files so we can resolve them in bulk
  private lazy val dxFiles: Vector[DxFile] = {
    val allInputs: Map[String, JsValue] = defaultsJs ++ inputsJs.values.flatten
    val fileJs = bundle.allCallables.values.toVector.flatMap { callable: Callable =>
      callable.inputVars.flatMap { param =>
        allInputs
          .get(s"${callable.name}.${param.name.decoded}")
          .map(jsv => extractDxFiles(jsv, param.dxType))
          .getOrElse(Vector.empty)
      }
    }
    val (dxFiles, dxPaths) = fileJs.foldLeft((Vector.empty[DxFile], Vector.empty[String])) {
      case ((files, paths), obj: JsObject) =>
        (files :+ DxFile.fromJson(dxApi, obj), paths)
      case ((files, paths), JsString(path)) =>
        (files, paths :+ path)
      case (_, other) =>
        throw new Exception(s"invalid file ${other}")
    }
    val resolvedPaths = dxApi
      .resolveDataObjectBulk(dxPaths, project)
      .map {
        case (key, dxFile: DxFile) => key -> dxFile
        case (_, dxobj) =>
          throw new Exception(s"Scanning the input file produced ${dxobj} which is not a file")
      }
    // lookup platform files in bulk
    dxApi.describeFilesBulk(dxFiles ++ resolvedPaths.values, validate = true)
  }

  private lazy val dxFileDescCache: DxFileDescCache = DxFileDescCache(dxFiles)
  lazy val fileResolver: FileSourceResolver = {
    val dxProtocol = DxFileAccessProtocol(dxFileCache = dxFileDescCache)
    baseFileResolver.replaceProtocol[DxFileAccessProtocol](dxProtocol)
  }
  private lazy val parameterLinkSerializer =
    ParameterLinkSerializer(fileResolver, dxApi, pathsAsObjects = complexPathValues)
  private lazy val parameterLinkDeserializer = ParameterLinkDeserializer(dxFileDescCache, dxApi)

  private def deserializationHandler(jsValue: JsValue, t: Type): Either[JsValue, Value] = {
    Left(translateJsInput(jsValue, t))
  }

  lazy val bundleWithDefaults: Bundle = {
    logger.trace(s"Embedding defaults into the IR")
    val defaultsExactlyOnce = ExactlyOnce("default", defaultsJs, logger)
    val allCallablesWithDefaults: Map[String, Callable] = bundle.allCallables.map {
      case (callableName, applet: Application) =>
        val inputsWithDefaults = applet.inputs.map { param =>
          val key = s"${applet.name}.${param.name.decoded}"
          defaultsExactlyOnce.get(key) match {
            case None => param
            case Some(default: JsValue) =>
              val irValue = parameterLinkDeserializer.deserializeInputWithType(
                  default,
                  param.dxType,
                  key,
                  Some(deserializationHandler)
              )
              param.copy(defaultValue = Some(irValue))
          }
        }
        callableName -> applet.copy(inputs = inputsWithDefaults)
      case (name, workflow: Workflow) =>
        val workflowWithDefaults =
          if (workflow.locked) {
            // locked workflow - we have workflow-level inputs
            val inputsWithDefaults = workflow.inputs.map {
              case (param, stageInput) =>
                val key = s"${workflow.name}.${param.name.decoded}"
                val stageInputWithDefault = defaultsExactlyOnce.get(key) match {
                  case None => stageInput
                  case Some(default: JsValue) =>
                    val irValue = parameterLinkDeserializer.deserializeInputWithType(
                        default,
                        param.dxType,
                        key,
                        Some(deserializationHandler)
                    )
                    StaticInput(irValue)
                }
                (param, stageInputWithDefault)
            }
            workflow.copy(inputs = inputsWithDefaults)
          } else {
            // Workflow is unlocked, we don't have workflow-level inputs.
            // Instead, set the defaults in the common stage.
            val stagesWithDefaults = workflow.stages.map { stage =>
              val callee: Callable = bundle.allCallables(stage.calleeName)
              logger.trace(s"addDefaultToStage ${stage.dxStage.id}, ${stage.description}")
              val prefix = if (stage.dxStage.id == s"stage-${Constants.CommonStage}") {
                workflow.name
              } else {
                s"${workflow.name}.${stage.description}"
              }
              val inputsWithDefaults = stage.inputs.zipWithIndex.map {
                case (stageInput, idx) =>
                  val param = callee.inputVars(idx)
                  val key = s"${prefix}.${param.name.decoded}"
                  defaultsExactlyOnce.get(key) match {
                    case None => stageInput
                    case Some(default: JsValue) =>
                      val irValue =
                        parameterLinkDeserializer.deserializeInputWithType(
                            default,
                            param.dxType,
                            key,
                            Some(deserializationHandler)
                        )
                      StaticInput(irValue)
                  }
              }
              stage.copy(inputs = inputsWithDefaults)
            }
            workflow.copy(stages = stagesWithDefaults)
          }
        // check that the stage order hasn't changed
        val allStageNames = workflow.stages.map(_.dxStage)
        val embedAllStageNames = workflowWithDefaults.stages.map(_.dxStage)
        assert(allStageNames == embedAllStageNames)
        name -> workflowWithDefaults
      case other =>
        throw new Exception(s"Unexpected callable ${other}")
    }
    val primaryCallableWithDefaults =
      bundle.primaryCallable.map(primary => allCallablesWithDefaults(primary.name))
    defaultsExactlyOnce.checkAllUsed()
    bundle.copy(primaryCallable = primaryCallableWithDefaults,
                allCallables = allCallablesWithDefaults)

  }

  private lazy val tasks: Vector[Application] = bundle.allCallables.collect {
    case (_, callable: Application) => callable
  }.toVector

  // CWL packed workflows always use 'main' as the main process name -
  // this enables CwlInputTranslator to replace it with the filename.
  // TODO: the right solution to this is to have the pre-processing
  //  script change #main to #<filename>
  protected val mainPrefix: Option[String] = None

  private def translateJsonInputs(fields: Map[String, JsValue]): Map[String, (Type, Value)] = {
    val fieldsExactlyOnce = ExactlyOnce("input", fields, logger)

    def checkAndBindCallableInputs(callable: Callable,
                                   inputPrefix: Option[String] = None,
                                   dxPrefix: Option[String] = None): Map[String, (Type, Value)] = {
      val lookupPrefix = inputPrefix.getOrElse(callable.name)
      callable.inputVars.flatMap { parameter =>
        // the key for looking up the parameter value
        val key = s"${lookupPrefix}.${parameter.name.decoded}"
        // the input name for the dx applet/workflow
        val fieldName = if (useManifests) {
          parameter.name.decoded
        } else {
          parameter.name.encoded
        }
        // the dx input name
        val dxName = dxPrefix.map(p => s"${p}.${fieldName}").getOrElse(fieldName)
        fieldsExactlyOnce.get(key) match {
          case None        => Map.empty
          case Some(value) =>
            // Do not assign the value to any later stages. We found the variable
            // declaration, the others are variable uses.
            logger.trace(s"checkAndBind, found: ${key} -> ${dxName}")
            val irValue = parameterLinkDeserializer.deserializeInputWithType(
                value,
                parameter.dxType,
                key,
                Some(deserializationHandler)
            )
            Map(dxName -> (parameter.dxType, irValue))
        }
      }.toMap
    }

    val inputs: Map[String, (Type, Value)] = bundle.primaryCallable match {
      // File with WDL tasks only, no workflows
      case None if tasks.isEmpty => Map.empty
      case None if tasks.size > 1 =>
        throw new Exception(s"cannot generate one input file for ${tasks.size} tasks")
      case None                            => checkAndBindCallableInputs(tasks.head)
      case Some(app: Application)          => checkAndBindCallableInputs(app, inputPrefix = mainPrefix)
      case Some(wf: Workflow) if wf.locked =>
        // locked workflow - user can only set workflow-level inputs
        checkAndBindCallableInputs(wf, inputPrefix = mainPrefix)
      case Some(wf: Workflow) if useManifests =>
        throw new Exception(s"cannot use manifests with unlocked workflow ${wf.name}")
      case Some(wf: Workflow) =>
        // unlocked workflow with at least one stage - inputs go into the common stage
        val commonStage = wf.stages.head.dxStage.id
        val commonInputs = checkAndBindCallableInputs(wf, dxPrefix = Some(commonStage))
        // filter out auxiliary stages
        val auxStages = Set(
            Constants.CommonStage,
            Constants.EvalStage,
            Constants.OutputStage,
            Constants.ReorgStage
        ).map(stg => s"stage-${stg}")
        val middleStages = wf.stages.filterNot(stg => auxStages.contains(stg.dxStage.id))
        // Inputs for top level calls
        val middleInputs = middleStages.flatMap { stage =>
          // Find the input definitions for the stage, by locating the callee
          val callee: Callable = bundle.allCallables.getOrElse(
              stage.calleeName,
              throw new Exception(s"callable ${stage.calleeName} is missing")
          )
          checkAndBindCallableInputs(callee,
                                     inputPrefix = Some(s"${wf.name}.${stage.description}"),
                                     dxPrefix = Some(stage.description))
        }
        commonInputs ++ middleInputs
      case other => throw new Exception(s"Unknown case ${other.getClass}")
    }
    if (!ignoreUnusedInputs) {
      fieldsExactlyOnce.checkAllUsed()
    }
    inputs
  }

  /**
    * Build a dx input file, based on the raw input file and the Bundle.
    * The general idea is to figure out the ancestry of each app(let)/call/workflow
    * input. This provides the fully-qualified-name (fqn) of each IR variable. Then
    * we check if the fqn is defined in the input file.
    * @return
    */
  lazy val translatedInputs: Map[Path, Map[String, (Type, Value)]] = {
    inputsJs.map {
      case (path, inputs) =>
        logger.trace(s"Translating input file ${path}")
        path -> translateJsonInputs(inputs)
    }
  }

  def writeTranslatedInputs(): Unit = {
    translatedInputs.foreach {
      case (path, inputs) =>
        val jsValues = if (useManifests) {
          val (types, values) = inputs.map {
            case (name, (t, v)) =>
              val dxName = dxNameFactory.fromDecodedName(name)
              (dxName -> t, dxName -> ValueSerde.serializeWithType(v, t))
          }.unzip
          val manifest = Manifest(values.toMap, Some(types.toMap))
          JsObject(Constants.InputManifest.encoded -> manifest.toJson())
        } else {
          JsObject(inputs.flatMap {
            case (name, (t, v)) =>
              parameterLinkSerializer.createFields(dxNameFactory.fromDecodedName(name), t, v).map {
                case (dxName, jsv) => dxName.encoded -> jsv
              }
          })
        }
        val fileName = FileUtils.replaceFileSuffix(path, ".dx.json")
        val dxInputFile = path.getParent match {
          case null   => Paths.get(fileName)
          case parent => parent.resolve(fileName)
        }
        logger.trace(s"Writing DNAnexus JSON input file ${dxInputFile}")
        JsUtils.jsToFile(jsValues, dxInputFile)
    }
  }

  def writeTranslatedInputManifest(): Unit = {
    translatedInputs.foreach {
      case (path, inputs) =>
        val fileName = FileUtils.replaceFileSuffix(path, ".manifest.json")
        val manifestFile = path.getParent match {
          case null   => Paths.get(fileName)
          case parent => parent.resolve(fileName)
        }
        val (types, values) = inputs.map {
          case (name, (t, v)) =>
            val dxName = dxNameFactory.fromDecodedName(name)
            (dxName -> t, dxName -> ValueSerde.serializeWithType(v, t))
        }.unzip
        val manifest = Manifest(values.toMap, Some(types.toMap))
        JsUtils.jsToFile(manifest.toJson(), manifestFile)
    }
  }
}
