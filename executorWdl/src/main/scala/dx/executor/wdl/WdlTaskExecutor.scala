package dx.executor.wdl

import dx.api.{DxPath, InstanceTypeRequest}
import dx.core.io.StreamFiles
import dx.core.ir.{DxName, Type, Value, ValueSerde}
import dx.core.languages.wdl.{
  IrToWdlValueBindings,
  Runtime,
  VersionSupport,
  WdlDxName,
  WdlOptions,
  WdlUtils
}
import dx.executor.{JobMeta, TaskExecutor}
import dx.util.{Bindings, DockerUtils, Logger}
import spray.json.JsObject
import wdlTools.eval.WdlValues._
import wdlTools.eval.{Eval, WdlValueBindings}
import wdlTools.exec.{TaskCommandFileGenerator, TaskInputOutput}
import wdlTools.types.WdlTypes._
import wdlTools.types.{TypedAbstractSyntax => TAT}

object WdlTaskExecutor {
  def create(
      jobMeta: JobMeta,
      streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
      checkInstanceType: Boolean
  ): WdlTaskExecutor = {
    val wdlOptions = jobMeta.parserOptions.map(WdlOptions.fromJson).getOrElse(WdlOptions.default)
    logger.warning(wdlOptions.toString())
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceString(jobMeta.sourceCode, wdlOptions, jobMeta.fileResolver)
    if (doc.workflow.isDefined) {
      throw new Exception("a workflow shouldn't be a member of this document")
    }
    val tasks = doc.elements.collect {
      case task: TAT.Task => task.name -> task
    }.toMap
    logger.warning(tasks.toString())
    if (tasks.isEmpty) {
      throw new Exception("no tasks in this WDL program")
    }
    if (tasks.size > 1) {
      throw new Exception("More than one task in this WDL program")
    }
    WdlTaskExecutor(tasks.values.head,
                    versionSupport,
                    typeAliases,
                    jobMeta,
                    streamFiles,
                    checkInstanceType)
  }
}

case class WdlTaskExecutor(task: TAT.Task,
                           versionSupport: VersionSupport,
                           typeAliases: Bindings[String, T_Struct],
                           jobMeta: JobMeta,
                           streamFiles: StreamFiles.StreamFiles,
                           checkInstanceType: Boolean)
    extends TaskExecutor(jobMeta, streamFiles, checkInstanceType) {

  private val fileResolver = jobMeta.fileResolver
  private val logger = jobMeta.logger
  private lazy val evaluator = Eval(
      jobMeta.workerPaths,
      Some(versionSupport.version),
      Vector.empty,
      jobMeta.fileResolver,
      Logger.Quiet
  )
  private lazy val taskIO = TaskInputOutput(task, logger)

  override val executorName: String = "dxExecutorWdl"

  override protected lazy val getSchemas: Map[String, Type.TSchema] = {
    typeAliases.toMap.view.mapValues(WdlUtils.toIRSchema).toMap
  }

  private lazy val inputTypes: Map[DxName, T] = {
    task.inputs.map(d => WdlDxName.fromSourceName(d.name) -> d.wdlType).toMap
  }

  private def wdlInputs: Map[DxName, V] = {
    // convert IR to WDL values; discard auxiliary fields
    val inputWdlValues: Map[DxName, V] = jobMeta.primaryInputs.map {
      case (dxName, value) =>
        dxName -> WdlUtils.fromIRValue(value, inputTypes(dxName), dxName.decoded)
    }
    // add default values for any missing inputs
    // Enable special handling for unset array values -
    // DNAnexus does not distinguish between null and empty for
    // array inputs, so we treat a null value for a non-optional
    // array that is allowed to be empty as the empty array.
    trace("Evaluating default values for inputs")
    logger.warning(inputWdlValues.toString())
    val wdlInputs: Map[DxName, V] = taskIO
      .inputsFromValues(inputWdlValues.map {
        case (dxName, v) => dxName.decoded -> v
      }, evaluator, ignoreDefaultEvalError = false, nullCollectionAsEmpty = true)
      .toMap
      .map {
        case (name, v) => WdlDxName.fromSourceName(name) -> v
      }
    if (logger.isVerbose) {
        val inputStr = task.inputs
          .map { inputDef =>
            val value = wdlInputs.get(WdlDxName.fromSourceName(inputDef.name))
            s"${inputDef.name} -> (${inputDef.wdlType}, ${value})"
          }
          .mkString("\n  ")
        logger.traceLimited(s"WDL inputs:\n  ${inputStr}")
    }
    wdlInputs
  }

  override protected def getInputsWithDefaults: Map[DxName, (Type, Value)] = {
    WdlUtils.toIR(wdlInputs.map {
      case (k, v) => k -> (inputTypes(k), v)
    })
  }

  private def evaluatePrivateVariables(inputs: Map[DxName, V]): Map[DxName, V] = {
    // evaluate the private variables using the inputs
    val init: Bindings[String, V] = WdlValueBindings(inputs.map {
      case (dxName, v) => dxName.decoded -> v
    })
    task.privateVariables
      .foldLeft(init) {
        case (env, TAT.PrivateVariable(name, wdlType, expr)) =>
          val wdlValue = evaluator.applyExprAndCoerce(expr, wdlType, env)
          env.add(name, wdlValue)
      }
      .toMap
      .map {
        case (name, v) => WdlDxName.fromSourceName(name) -> v
      }
  }

  lazy val (runtimeOverrides, hintOverrides) = {
    val (runtimeOverridesJs, hintOverridesJs) = jobMeta.jsOverrides match {
      case Some(JsObject(fields)) if fields.contains("runtime") || fields.contains("hints") =>
        (fields.get("runtime").map(_.asJsObject.fields),
         fields.get("hints").map(_.asJsObject.fields))
      case Some(JsObject(fields)) => (Some(fields), None)
      case None                   => (None, None)
      case Some(other) =>
        throw new Exception(s"invalid overrides ${other}")
    }
    (runtimeOverridesJs.map(o => IrToWdlValueBindings(ValueSerde.deserializeMap(o))),
     hintOverridesJs.map(o => IrToWdlValueBindings(ValueSerde.deserializeMap(o))))
  }

  private def createRuntime(env: Map[DxName, V]): Runtime = {
    Runtime(
        versionSupport.version,
        task.runtime,
        task.hints,
        evaluator,
        runtimeOverrides,
        hintOverrides,
        Some(IrToWdlValueBindings(jobMeta.defaultRuntimeAttrs)),
        Some(WdlValueBindings(env.map {
          case (dxName, v) => dxName.decoded -> v
        }))
    )
  }

  override protected def getInstanceTypeRequest(
      inputs: Map[DxName, (Type, Value)]
  ): InstanceTypeRequest = {
    val wdlInputs = WdlUtils.fromIR(inputs, typeAliases.toMap).map {
      case (name, (_, value)) => name -> value
    }
    val env = evaluatePrivateVariables(wdlInputs)
    val runtime = createRuntime(env)
    runtime.parseInstanceType
  }

  private lazy val hints =
    HintResolver(versionSupport.version, task.parameterMeta, task.hints, hintOverrides)

  /**
    * Should we try to stream the file(s) associated with the given input parameter?
    * This can be set at the parameter level (in parameters_meta or hints.inputs) or
    * at the global level (at the hints top level).
    */
  override protected def streamFileForInput(parameterName: DxName): Boolean = {
    hints.isLocalizationOptional(parameterName.decoded)
  }

  override protected def writeCommandScript(
      localizedInputs: Map[DxName, (Type, Value)],
      localizedDependencies: Option[Map[String, (Type, Value)]]
  ): (Map[DxName, (Type, Value)], Boolean, Option[Set[Int]]) = {
    val inputs = WdlUtils.fromIR(localizedInputs, typeAliases.toMap)
    val inputValues = inputs.map {
      case (name, (_, v)) => name -> v
    }
    val inputsWithPrivateVars = evaluatePrivateVariables(inputValues)
    val inputAndPrivateVarTypes = inputTypes ++ task.privateVariables
      .map(d => WdlDxName.fromSourceName(d.name) -> d.wdlType)
      .toMap
    val updatedInputs = WdlUtils.toIR(inputsWithPrivateVars.map {
      case (name, value) => name -> (inputAndPrivateVarTypes(name), value)
    })
    evaluator.applyCommand(task.command, WdlValueBindings(inputsWithPrivateVars.map {
      case (dxName, v) => dxName.decoded -> v
    })) match {
      case command if command.trim.isEmpty => (updatedInputs, false, None)
      case command =>
        val generator = TaskCommandFileGenerator(logger)
        val runtime = createRuntime(inputsWithPrivateVars)
        val dockerUtils = DockerUtils(fileResolver, logger)
        val container = runtime.container match {
          case Vector() => None
          case Vector(image) =>
            val resolvedImage = dockerUtils.getImage(image)
            Some(resolvedImage, jobMeta.workerPaths)
          case images =>
            // we prefer a dx:// url, false comes before true in sort order
            val resolvedImage =
              dockerUtils.getImage(images.sortBy(!_.startsWith(DxPath.DxUriPrefix)))
            Some(resolvedImage, jobMeta.workerPaths)
        }
        generator.apply(Some(command), jobMeta.workerPaths, container)
        (updatedInputs, true, runtime.returnCodes)
    }
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[DxName, (Type, Value)]
  ): (Map[DxName, (Type, Value)], Map[DxName, (Set[String], Map[String, String])]) = {
    val outputTypes: Map[String, T] = task.outputs.map(d => d.name -> d.wdlType).toMap
    // Evaluate the output parameters in dependency order.
    val localizedOutputs = WdlUtils.toIR(
        taskIO
          .evaluateOutputs(
              evaluator,
              WdlValueBindings(
                  WdlUtils.fromIR(localizedInputs, typeAliases.toMap).map {
                    case (dxName, (_, v)) => dxName.decoded -> v
                  }
              )
          )
          .toMap
          .map {
            case (name, value) =>
              WdlDxName.fromSourceName(name) -> (outputTypes(name), value)
          }
    )
    val tagsAndProperties = localizedOutputs.keys.map { name =>
      val (tags, properties) = hints.getOutput(name.decoded) match {
        case Some(V_Object(fields)) =>
          val tags = fields.get("tags") match {
            case Some(V_Array(tags)) =>
              tags.map {
                case V_String(tag) => tag
                case other         => throw new Exception(s"invalid tag ${other}")
              }.toSet
            case _ => Set.empty[String]
          }
          val properties = fields.get("properties") match {
            case Some(V_Object(properties)) =>
              properties.map {
                case (key, V_String(value)) => key -> value
                case other                  => throw new Exception(s"invalid property ${other}")
              }
            case _ => Map.empty[String, String]
          }
          (tags, properties)
        case _ => (Set.empty[String], Map.empty[String, String])
      }
      name -> (tags, properties)
    }.toMap
    (localizedOutputs, tagsAndProperties)
  }

  override protected lazy val outputTypes: Map[DxName, Type] = {
    task.outputs.map { outputDef: TAT.OutputParameter =>
      WdlDxName.fromSourceName(outputDef.name) -> WdlUtils.toIRType(outputDef.wdlType)
    }.toMap
  }
}
