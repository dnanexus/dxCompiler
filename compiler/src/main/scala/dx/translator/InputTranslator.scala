package dx.translator

import java.nio.file.{Path, Paths}
import dx.api.{DxApi, DxFile, DxFileDescCache, DxPath, DxProject}
import dx.core.Constants
import dx.core.ir.Type._
import dx.core.ir.{
  Application,
  Bundle,
  Callable,
  DxName,
  DxNameFactory,
  Manifest,
  ParameterLinkDeserializer,
  ParameterLinkSerializer,
  StageInputStatic,
  Type,
  Value,
  ValueSerde,
  Workflow
}
import dx.util.{FileSourceResolver, FileUtils, JsUtils, Logger}
import dx.util.CollectionUtils.IterableOnceExtensions
import dx.util.protocols.DxFileAccessProtocol
import dx.yaml._
import spray.json._

import scala.util.Try

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

object InputTranslator {
  val ManifestSuffix = s".${Constants.InputManifest.decoded}"

  def loadJsonFileWithComments(path: Path): Map[String, JsValue] = {
    // skip comment lines, which start with ##
    JsUtils
      .getFields(JsUtils.jsFromFile(path))
      .view
      .filterKeys(!_.startsWith("##"))
      .toMap
  }

  private def yamlToJson(yamlValue: YamlValue): JsValue = {
    yamlValue match {
      case YamlNull        => JsNull
      case YamlString(s)   => JsString(s)
      case YamlBoolean(b)  => JsBoolean(b)
      case YamlNumber(n)   => JsNumber(n)
      case YamlNaN         => JsNumber(Double.NaN)
      case YamlPositiveInf => JsNumber(Double.PositiveInfinity)
      case YamlNegativeInf => JsNumber(Double.NegativeInfinity)
      case YamlArray(a)    => JsArray(a.map(yamlToJson))
      case YamlSet(s)      => JsArray(s.map(yamlToJson).toVector)
      case YamlObject(obj) =>
        JsObject(obj.map {
          case (YamlString(key), value) => key -> yamlToJson(value)
          case (key, value)             => key.toString -> yamlToJson(value)
        })
      case other =>
        throw new Exception(s"unexpected value ${other}")
    }
  }

  // load YAML and convert to JSON
  def loadYamlFile(path: Path): Map[String, JsValue] = {
    yamlToJson(FileUtils.readFileContent(path).parseYaml) match {
      case JsObject(fields) => fields
      case other =>
        throw new Exception(s"expected object not ${other}")
    }
  }

  def loadInputFile(path: Path): Map[String, JsValue] = {
    val filename = path.getFileName.toString
    if (filename.endsWith(".json")) {
      loadJsonFileWithComments(path)
    } else if (filename.endsWith(".yaml") || filename.endsWith(".yml")) {
      loadYamlFile(path)
    } else {
      Try(loadJsonFileWithComments(path)).orElse(Try(loadYamlFile(path))).get
    }
  }
}

abstract class InputTranslator(bundle: Bundle,
                               inputPaths: Vector[Path],
                               defaultPaths: Option[Path],
                               project: DxProject,
                               useManifests: Boolean,
                               complexPathValues: Boolean,
                               ignoreUnusedInputs: Boolean,
                               dxNameFactory: DxNameFactory,
                               baseFileResolver: FileSourceResolver = FileSourceResolver.get,
                               dxApi: DxApi = DxApi.get,
                               logger: Logger = Logger.get) {

  /**
    * Splits raw inputs into two maps: 1) parameter values and 2) overrides.
    * The overrides value is language-specific and is serialized as-is.
    */
  protected def splitInputs(
      rawInputs: Map[String, JsValue]
  ): (Map[String, JsValue], Map[String, JsObject])

  // All of the inputs and overrides in JSON format. Keys are of the form
  // <prefix>.<decoded_name>, where prefix is the task or workflow name and
  // decoded_name is the parameter name as it appears in the source.
  private lazy val (allRawInputs: Map[Path, Map[String, JsValue]],
                    allRawOverrides: Map[Path, Map[String, JsObject]]) = {
    val (rawInputs, overrides) = inputPaths.map { path =>
      val jsValues = InputTranslator.loadInputFile(path)
      val (inputs, overridesJs) = splitInputs(jsValues)
      val inputsJs =
        if (inputs.size == 1 && inputs.keys.head.endsWith(InputTranslator.ManifestSuffix)) {
          val (key, fields) = inputs.head
          val prefix = key.dropRight(InputTranslator.ManifestSuffix.length)
          val manifest = Manifest.parse(fields, dxNameFactory)
          manifest.jsValues.map {
            case (name, value) => s"${prefix}.${name}" -> value
          }
        } else {
          inputs
        }
      (path -> inputsJs, path -> overridesJs)
    }.unzip
    (rawInputs.toMap, overrides.toMap)
  }

  /**
    * Converts a language-specific JSON value to one that can be deserialized to an IR Value.
    * This is the mechanism to support compile-time inputs/defaults that do not map implicitly
    * to an IR type.
    */
  protected def convertRawInput(rawInput: JsValue, t: Type): JsValue = rawInput

  /**
    * Extracts Dx files from a JSON input.
    * @param rawInput the JSON value
    * @param t the parameter type
    * @return a Vector of JSON values that describe Dx files
    */
  private def extractDxFiles(rawInput: JsValue, t: Type): Vector[JsValue] = {
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

    val updatedValue = convertRawInput(rawInput, t)
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

  private lazy val defaults =
    defaultPaths.map(InputTranslator.loadJsonFileWithComments).getOrElse(Map.empty)

  // scan all inputs and defaults and extract files so we can resolve them in bulk
  private lazy val dxFiles: Vector[DxFile] = {
    val allInputs: Vector[Map[String, JsValue]] = Vector(defaults) ++ allRawInputs.values.toVector
    val fileJs = bundle.allCallables.values.toVector.flatMap { callable: Callable =>
      callable.inputVars.flatMap { param =>
        val key = s"${callable.name}.${param.name.decoded}"
        allInputs.flatMap(m =>
          m.get(key)
            .map(jsv => extractDxFiles(jsv, param.dxType))
            .getOrElse(Vector.empty)
        )
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
    Left(convertRawInput(jsValue, t))
  }

  lazy val bundleWithDefaults: Bundle = {
    if (defaults.nonEmpty) {
      logger.trace(s"Embedding defaults into the IR")
      val defaultsExactlyOnce = ExactlyOnce("default", defaults, logger)
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
                  logger.warning("inputtranslator: bundlewithDefaults: irvalue HERE!")
                  val stageInputWithDefault = defaultsExactlyOnce.get(key) match {
                    case None => stageInput
                    case Some(default: JsValue) =>
                      val irValue = parameterLinkDeserializer.deserializeInputWithType(
                          default,
                          param.dxType,
                          key,
                          Some(deserializationHandler)
                      )
                      logger.warning("irValue")
                      logger.warning(irValue.toString())
                      StageInputStatic(irValue)
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
                    logger.warning("inputtranslator2: bundlewithDefaults2: irvalue HERE!")
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
                        logger.warning("irValue2")
                        logger.warning(irValue.toString())
                        StageInputStatic(irValue)
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

    } else {
      bundle
    }
  }

  private lazy val applications: Vector[Application] = bundle.allCallables.collect {
    case (_, callable: Application) => callable
  }.toVector

  /**
    * The prefix to add to the field name when looking up the field
    * in the inputs.
    */
  protected val mainPrefix: Option[String] = None

  /**
    * Whether it is optional for the input to a locked workflow to
    * be prefixed with the workflow name.
    */
  protected val lockedWorkflowPrefixOptional: Boolean = false

  // type aliases for translated input and override info
  private type DxInput = (Option[String], DxName, Type, Value)
  private type DxOverride = (Option[String], JsObject)

  /**
    * Translates language-specific inputs and overrides to dx format. The general
    * idea is to figure out the ancestry of each app(let)/call/workflow input.
    * This provides the fully-qualified-name (fqn) of each IR variable. Then we
    * check if the fqn is defined in the input file.
    * @param rawInputs inputs to translate
    * @return Vector of (prefix, name, type, value), where prefix is optional.
    */
  private def translate(
      rawInputs: Map[String, JsValue],
      overrides: Option[Map[String, JsObject]]
  ): (Vector[(Option[String], DxName, Type, Value)], Vector[DxOverride]) = {
    val fieldsExactlyOnce = ExactlyOnce("input", rawInputs, logger)

    def bind(
        callable: Callable,
        inputPrefix: Option[String] = None,
        inputPrefixOptional: Boolean = false,
        dxPrefix: Option[String] = None
    ): (Vector[DxInput], Vector[DxOverride]) = {
      val lookupPrefix = inputPrefix.getOrElse(callable.name)
      val boundInputs = callable.inputVars.flatMap { parameter =>
        // the key(s) for looking up the parameter value
        val keys = Vector(
            Some(s"${lookupPrefix}.${parameter.name.decoded}"),
            Option.when(inputPrefixOptional)(parameter.name.decoded)
        ).flatten
        // the dx input name
        keys
          .collectFirstDefined(key => fieldsExactlyOnce.get(key).map((key, _)))
          .map {
            case (key, value) =>
              // Do not assign the value to any later stages. We found the variable
              // declaration, the others are variable uses.
              logger.trace(s"checkAndBind, found: ${key} -> ${parameter.name}")
              val irValue = parameterLinkDeserializer.deserializeInputWithType(
                  value,
                  parameter.dxType,
                  key,
                  Some(deserializationHandler)
              )
              (dxPrefix, parameter.name, parameter.dxType, irValue)
          }
          .orElse {
            if (!Type.isOptional(parameter.dxType) && parameter.defaultValue.isEmpty) {
              logger.warning(s"missing input for non-optional parameter ${parameter.name}")
            }
            None
          }
      }
      val boundOverrides = overrides
        .flatMap(_.get(lookupPrefix).map { prefixOverrides =>
          (dxPrefix, JsObject("___" -> prefixOverrides))
        })
        .toVector
      (boundInputs, boundOverrides)
    }

    val (dxInputs, dxOverrides) = {
      bundle.primaryCallable match {
        case None if applications.isEmpty => (Vector.empty[DxInput], Vector.empty[DxOverride])
        case None if applications.size > 1 =>
          throw new Exception(s"cannot generate one input file for ${applications.size} tasks")
        case None =>
          bind(applications.head, inputPrefix = mainPrefix, inputPrefixOptional = true)
        case Some(app: Application) =>
          bind(app, inputPrefix = mainPrefix, inputPrefixOptional = true)
        case Some(wf: Workflow) if wf.locked =>
          // locked workflow - user can only set workflow-level inputs
          bind(wf, inputPrefix = mainPrefix, inputPrefixOptional = lockedWorkflowPrefixOptional)
        case Some(wf: Workflow) if useManifests =>
          throw new Exception(s"cannot use manifests with unlocked workflow ${wf.name}")
        case Some(wf: Workflow) =>
          // unlocked workflow with at least one stage - inputs go into the common stage
          val commonStage = wf.stages.head.dxStage.id
          val (commonInputs, commonOverrides) = bind(wf, dxPrefix = Some(commonStage))
          // filter out auxiliary stages
          val auxStages = Set(
              Constants.CommonStage,
              Constants.EvalStage,
              Constants.OutputStage,
              Constants.ReorgStage
          ).map(stg => s"stage-${stg}")
          val middleStages = wf.stages.filterNot(stg => auxStages.contains(stg.dxStage.id))
          // Inputs for top level calls
          val (middleInputs, middleOverrides) = middleStages.map { stage =>
            // Find the input definitions for the stage, by locating the callee
            val callee: Callable = bundle.allCallables.getOrElse(
                stage.calleeName,
                throw new Exception(s"callable ${stage.calleeName} is missing")
            )
            bind(callee,
                 inputPrefix = Some(s"${wf.name}.${stage.description}"),
                 dxPrefix = Some(stage.description))
          }.unzip
          (commonInputs ++ middleInputs.flatten, commonOverrides ++ middleOverrides.flatten)
        case other => throw new Exception(s"Unknown case ${other.getClass}")
      }
    }
    if (!ignoreUnusedInputs) {
      fieldsExactlyOnce.checkAllUsed()
    }
    (dxInputs, dxOverrides)
  }

  def writeTranslatedInputs(): Unit = {
    allRawInputs.foreach {
      case (path, inp) =>
        logger.trace(s"Translating input file ${path}")
        val (dxInputs, dxOverrides) = translate(inp, allRawOverrides.get(path))
        val jsValues = if (useManifests) {
          val (types, values) = dxInputs.map {
            case (_, dxName, t, v) =>
              (dxName -> t, dxName -> ValueSerde.serializeWithType(v, t))
          }.unzip
          if (dxOverrides.size > 1) {
            throw new Exception(
                s"manifest cannot have overrides for multiple executables: ${dxOverrides}"
            )
          }
          val (overrideTypes, overrideValues) = dxOverrides.headOption
            .map {
              case (_, overridesJs) =>
                (Map(Constants.Overrides -> THash), Map(Constants.Overrides -> overridesJs))
            }
            .getOrElse((Map.empty[DxName, Type], Map.empty[DxName, JsValue]))
          val manifest =
            Manifest(values.toMap ++ overrideValues, Some(types.toMap ++ overrideTypes))
          JsObject(Constants.InputManifest.encoded -> manifest.toJson())
        } else {
          val translatedInputs = dxInputs.flatMap {
            case (prefix, dxName, t, v) =>
              parameterLinkSerializer.createFields(dxName, t, v).map {
                case (dxName, jsv) =>
                  val inputName =
                    prefix.map(p => s"${p}.${dxName.encoded}").getOrElse(dxName.encoded)
                  inputName -> jsv
              }
          }.toMap
          val translatedOverrides = dxOverrides.map {
            case (prefix, overridesJs) =>
              prefix
                .map(p => s"${p}.${Constants.Overrides.encoded}")
                .getOrElse(Constants.Overrides.encoded) -> overridesJs
          }.toMap
          JsObject(translatedInputs ++ translatedOverrides)
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
    allRawInputs.foreach {
      case (path, inp) =>
        logger.trace(s"Translating input file ${path}")
        val (dxInputs, dxOverrides) = translate(inp, allRawOverrides.get(path))
        val fileName = FileUtils.replaceFileSuffix(path, ".manifest.json")
        val manifestFile = path.getParent match {
          case null   => Paths.get(fileName)
          case parent => parent.resolve(fileName)
        }
        val (types, values) = dxInputs.map {
          case (_, dxName, t, v) =>
            (dxName -> t, dxName -> ValueSerde.serializeWithType(v, t))
        }.unzip
        if (dxOverrides.size > 1) {
          throw new Exception(
              s"manifest cannot have overrides for multiple executables: ${dxOverrides}"
          )
        }
        val (overrideTypes, overrideValues) = dxOverrides.headOption
          .map {
            case (_, overridesJs) =>
              (Map(Constants.Overrides -> THash), Map(Constants.Overrides -> overridesJs))
          }
          .getOrElse((Map.empty[DxName, Type], Map.empty[DxName, JsValue]))
        val manifest = Manifest(values.toMap ++ overrideValues, Some(types.toMap ++ overrideTypes))
        JsUtils.jsToFile(manifest.toJson(), manifestFile)
    }
  }
}
