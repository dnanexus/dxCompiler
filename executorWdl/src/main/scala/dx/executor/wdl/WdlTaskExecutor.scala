package dx.executor.wdl

import dx.api.{DxPath, InstanceTypeRequest}
import dx.core.io.StreamFiles
import dx.core.ir.{Type, Value}
import dx.core.languages.wdl.{IrToWdlValueBindings, Runtime, VersionSupport, WdlOptions, WdlUtils}
import dx.executor.{JobMeta, TaskExecutor}
import dx.util.{Bindings, DockerUtils, Logger}
import wdlTools.eval.WdlValues._
import wdlTools.eval.{Eval, WdlValueBindings}
import wdlTools.exec.{TaskCommandFileGenerator, TaskInputOutput}
import wdlTools.types.WdlTypes._
import wdlTools.types.{TypedAbstractSyntax => TAT}

object WdlTaskExecutor {
  def create(
      jobMeta: JobMeta,
      streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
      waitOnUpload: Boolean
  ): WdlTaskExecutor = {
    val wdlOptions = jobMeta.parserOptions.map(WdlOptions.fromJson).getOrElse(WdlOptions.default)
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceString(jobMeta.sourceCode, wdlOptions, jobMeta.fileResolver)
    if (doc.workflow.isDefined) {
      throw new Exception("a workflow shouldn't be a member of this document")
    }
    val tasks = doc.elements.collect {
      case task: TAT.Task => task.name -> task
    }.toMap
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
                    waitOnUpload)
  }
}

case class WdlTaskExecutor(task: TAT.Task,
                           versionSupport: VersionSupport,
                           typeAliases: Bindings[String, T_Struct],
                           jobMeta: JobMeta,
                           streamFiles: StreamFiles.StreamFiles,
                           waitOnUpload: Boolean)
    extends TaskExecutor(jobMeta, streamFiles, waitOnUpload) {

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

  private lazy val inputTypes: Map[String, T] = {
    task.inputs.map(d => d.name -> d.wdlType).toMap
  }

  private def wdlInputs: Map[String, V] = {
    // convert IR to WDL values; discard auxiliary fields
    val inputWdlValues: Map[String, V] = jobMeta.primaryInputs.map {
      case (name, value) =>
        name -> WdlUtils.fromIRValue(value, inputTypes(name), name)
    }
    // add default values for any missing inputs
    // Enable special handling for unset array values -
    // DNAnexus does not distinguish between null and empty for
    // array inputs, so we treat a null value for a non-optional
    // array that is allowed to be empty as the empty array.
    trace("Evaluating default values for inputs")
    val wdlInputs = taskIO
      .inputsFromValues(inputWdlValues,
                        evaluator,
                        ignoreDefaultEvalError = false,
                        nullCollectionAsEmpty = true)
      .toMap
    if (logger.isVerbose) {
      if (logger.isVerbose) {
        val inputStr = task.inputs
          .map { inputDef =>
            s"${inputDef.name} -> (${inputDef.wdlType}, ${wdlInputs.get(inputDef.name)})"
          }
          .mkString("\n  ")
        logger.traceLimited(s"WDL inputs:\n  ${inputStr}")
      }
    }
    wdlInputs
  }

  override protected def getInputsWithDefaults: Map[String, (Type, Value)] = {
    WdlUtils.toIR(wdlInputs.map {
      case (k, v) => k -> (inputTypes(k), v)
    })
  }

  private def evaluatePrivateVariables(inputs: Map[String, V]): Map[String, V] = {
    // evaluate the private variables using the inputs
    val env: Map[String, V] =
      task.privateVariables.foldLeft(inputs) {
        case (env, TAT.PrivateVariable(name, wdlType, expr)) =>
          val wdlValue =
            evaluator.applyExprAndCoerce(expr, wdlType, WdlValueBindings(env))
          env + (name -> wdlValue)
      }
    env
  }

  private def createRuntime(env: Map[String, V]): Runtime = {
    Runtime(
        versionSupport.version,
        task.runtime,
        task.hints,
        evaluator,
        Some(IrToWdlValueBindings(jobMeta.defaultRuntimeAttrs)),
        Some(WdlValueBindings(env))
    )
  }

  private def getRequiredInstanceTypeRequest(
      inputs: Map[String, V] = wdlInputs
  ): InstanceTypeRequest = {
    val env = evaluatePrivateVariables(inputs)
    val runtime = createRuntime(env)
    runtime.parseInstanceType
  }

  override protected lazy val getInstanceTypeRequest: InstanceTypeRequest =
    getRequiredInstanceTypeRequest()

  private lazy val hints =
    HintResolver(versionSupport.version, task.parameterMeta, task.hints)

  /**
    * Should we try to stream the file(s) associated with the given input parameter?
    * This can be set at the parameter level (in parameters_meta or hints.inputs) or
    * at the global level (at the hints top level).
    */
  override protected def streamFileForInput(parameterName: String): Boolean = {
    hints.isLocalizationOptional(parameterName)
  }

  override protected def writeCommandScript(
      localizedInputs: Map[String, (Type, Value)]
  ): (Map[String, (Type, Value)], Boolean, Option[Set[Int]], Set[Int]) = {
    val inputs = WdlUtils.fromIR(localizedInputs, typeAliases.toMap)
    val inputValues = inputs.map {
      case (name, (_, v)) => name -> v
    }
    val inputsWithPrivateVars = evaluatePrivateVariables(inputValues)
    val inputAndPrivateVarTypes = inputTypes ++ task.privateVariables
      .map(d => d.name -> d.wdlType)
      .toMap
    val updatedInputs = WdlUtils.toIR(inputsWithPrivateVars.map {
      case (name, value) => name -> (inputAndPrivateVarTypes(name), value)
    })
    evaluator.applyCommand(task.command, WdlValueBindings(inputsWithPrivateVars)) match {
      case command if command.trim.isEmpty =>
        (updatedInputs, false, None, Set.empty[Int])
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
        (updatedInputs, true, runtime.returnCodes, Set.empty[Int])
    }
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): (Map[String, (Type, Value)], Map[String, (Set[String], Map[String, String])]) = {
    val outputTypes: Map[String, T] = task.outputs.map(d => d.name -> d.wdlType).toMap
    // Evaluate the output parameters in dependency order.
    val localizedOutputs = WdlUtils.toIR(
        taskIO
          .evaluateOutputs(
              evaluator,
              WdlValueBindings(
                  WdlUtils.fromIR(localizedInputs, typeAliases.toMap).view.mapValues(_._2).toMap
              )
          )
          .toMap
          .map {
            case (name, value) => name -> (outputTypes(name), value)
          }
    )
    val tagsAndProperties = localizedOutputs.keys.map { name =>
      val (tags, properties) = hints.getOutput(name) match {
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

  override protected lazy val outputTypes: Map[String, Type] = {
    task.outputs.map { outputDef: TAT.OutputParameter =>
      outputDef.name -> WdlUtils.toIRType(outputDef.wdlType)
    }.toMap
  }
}
