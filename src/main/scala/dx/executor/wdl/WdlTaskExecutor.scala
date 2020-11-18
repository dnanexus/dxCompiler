package dx.executor.wdl

import dx.api.DxPath
import dx.core.{getVersion => getDxWdlVersion}
import dx.core.io.StreamFiles
import dx.core.ir.{Type, Value}
import dx.core.languages.wdl.{DxMetaHints, Runtime, VersionSupport, WdlUtils}
import dx.executor.{FileUploader, JobMeta, SerialFileUploader, TaskExecutor}
import dx.translator.wdl.IrToWdlValueBindings
import dx.util.{Bindings, Logger, TraceLevel}
import wdlTools.eval.WdlValues._
import wdlTools.eval.{Eval, Hints, Meta, WdlValueBindings}
import wdlTools.exec.{DockerUtils, TaskCommandFileGenerator, TaskInputOutput}
import wdlTools.syntax.SourceLocation
import wdlTools.types.TypeCheckingRegime.TypeCheckingRegime
import wdlTools.types.WdlTypes._
import wdlTools.types.{TypeCheckingRegime, TypedAbstractSyntax => TAT}

object WdlTaskExecutor {
  val MaxDisambiguationDirs: Int = 5000

  def create(
      jobMeta: JobMeta,
      fileUploader: FileUploader = SerialFileUploader(),
      streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
      regime: TypeCheckingRegime = TypeCheckingRegime.Moderate
  ): WdlTaskExecutor = {
    val (doc, typeAliases, versionSupport) =
      VersionSupport.fromSourceString(jobMeta.sourceCode, jobMeta.fileResolver, regime)
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
                    fileUploader,
                    streamFiles)
  }
}

case class WdlTaskExecutor(task: TAT.Task,
                           versionSupport: VersionSupport,
                           typeAliases: Bindings[String, T_Struct],
                           jobMeta: JobMeta,
                           fileUploader: FileUploader,
                           streamFiles: StreamFiles.StreamFiles)
    extends TaskExecutor(jobMeta, fileUploader, streamFiles) {

  private val fileResolver = jobMeta.fileResolver

  private val logger = jobMeta.logger

  private lazy val evaluator = Eval(
      jobMeta.workerPaths,
      Some(versionSupport.version),
      jobMeta.fileResolver,
      Logger.Quiet
  )
  private lazy val taskIO = TaskInputOutput(task, logger)

  override val executorName: String = "dxExecutorWdl"

  override def getVersion: String = getDxWdlVersion

  override protected lazy val getSchemas: Map[String, Type.TSchema] = {
    typeAliases.toMap.view.mapValues(WdlUtils.toIRSchema).toMap
  }

  private lazy val inputTypes: Map[String, T] =
    task.inputs.map(d => d.name -> d.wdlType).toMap

  private def wdlInputs: Map[String, V] = {
    // convert IR to WDL values; discard auxiliary fields
    val inputWdlValues: Map[String, V] = jobMeta.primaryInputs.map {
      case (name, value) =>
        name -> WdlUtils.fromIRValue(value, inputTypes(name), name)
    }
    // add default values for any missing inputs
    taskIO.inputsFromValues(inputWdlValues, evaluator, strict = true).toMap
  }

  private def printInputs(inputs: Map[String, V]): Unit = {
    if (logger.isVerbose) {
      val inputStr = task.inputs
        .map { inputDef =>
          s"${inputDef.name} -> (${inputDef.wdlType}, ${inputs.get(inputDef.name)})"
        }
        .mkString("\n")
      logger.traceLimited(s"inputs: ${inputStr}")
    }
  }

  override protected def getInputsWithDefaults: Map[String, (Type, Value)] = {
    val inputs = wdlInputs
    printInputs(inputs)
    WdlUtils.toIR(inputs.map {
      case (k, v) => k -> (inputTypes(k), v)
    })
  }

  private def evaluatePrivateVariables(inputs: Map[String, V]): Map[String, V] = {
    // evaluate the private variables using the inputs
    val env: Map[String, V] =
      task.privateVariables.foldLeft(inputs) {
        case (env, TAT.PrivateVariable(name, wdlType, expr, _)) =>
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

  private def getRequiredInstanceType(inputs: Map[String, V] = wdlInputs): String = {
    logger.traceLimited("calcInstanceType", minLevel = TraceLevel.VVerbose)
    printInputs(inputs)
    val env = evaluatePrivateVariables(inputs)
    val runtime = createRuntime(env)
    val request = runtime.parseInstanceType
    logger.traceLimited(s"calcInstanceType $request")
    jobMeta.instanceTypeDb.apply(request).name
  }

  override protected lazy val getRequiredInstanceType: String = getRequiredInstanceType()

  private lazy val parameterMeta = Meta.create(versionSupport.version, task.parameterMeta)

  /**
    * Should we try to stream the file(s) associated with the given input parameter?
    */
  override protected def streamFileForInput(parameterName: String): Boolean = {
    parameterMeta.get(parameterName) match {
      case Some(V_String(DxMetaHints.ParameterMetaStream)) =>
        true
      case Some(V_Object(fields)) =>
        // This enables the stream annotation in the object form of metadata value, e.g.
        // bam_file : {
        //   stream : true
        // }
        // We also support two aliases, dx_stream and localizationOptional
        fields.view
          .filterKeys(
              Set(DxMetaHints.ParameterMetaStream,
                  DxMetaHints.ParameterHintStream,
                  Hints.LocalizationOptionalKey)
          )
          .values
          .exists {
            case V_Boolean(b) => b
            case _            => false
          }
      case _ => false
    }
  }

  override protected def writeCommandScript(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    val inputs = WdlUtils.fromIR(localizedInputs, typeAliases.toMap)
    val inputValues = inputs.map {
      case (name, (_, v)) => name -> v
    }
    printInputs(inputValues)
    val inputsWithPrivateVars = evaluatePrivateVariables(inputValues)
    // TODO: there may be private variables that reference files created by the
    //  command, or functions that depend on the execution of the command
    //  (e.g. stdout()). Split the private vars into those that need to be
    //  evaluated before vs after the command, and only evaluate the former here.
    val ctx = WdlValueBindings(inputsWithPrivateVars)
    val command = evaluator.applyCommand(task.command, ctx) match {
      case s if s.trim.isEmpty => None
      case s                   => Some(s)
    }
    val generator = TaskCommandFileGenerator(logger)
    val runtime = createRuntime(inputsWithPrivateVars)
    val dockerUtils = DockerUtils(fileResolver, logger)
    val container = runtime.container match {
      case Vector() => None
      case Vector(image) =>
        val resolvedImage = dockerUtils.getImage(image, SourceLocation.empty)
        Some(resolvedImage, jobMeta.workerPaths)
      case v =>
        // we prefer a dx:// url
        val (dxUrls, imageNames) = v.partition(_.startsWith(DxPath.DxUriPrefix))
        val resolvedImage = dockerUtils.getImage(dxUrls ++ imageNames, SourceLocation.empty)
        Some(resolvedImage, jobMeta.workerPaths)
    }
    generator.apply(command, jobMeta.workerPaths, container)
    val inputAndPrivateVarTypes = inputTypes ++ task.privateVariables
      .map(d => d.name -> d.wdlType)
      .toMap
    WdlUtils.toIR(ctx.bindings.map {
      case (name, value) => name -> (inputAndPrivateVarTypes(name), value)
    })
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    val outputTypes: Map[String, T] = task.outputs.map(d => d.name -> d.wdlType).toMap
    // Evaluate the output parameters in dependency order.
    // These will include output files without canonicalized paths, which is why we need
    // the following complex logic to match up local outputs to remote URIs.
    // TODO: evaluate any private variables that depend on the command
    val localizedOutputs = taskIO
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
    WdlUtils.toIR(localizedOutputs)
  }

  override protected lazy val outputTypes: Map[String, Type] = {
    task.outputs.map { outputDef: TAT.OutputParameter =>
      outputDef.name -> WdlUtils.toIRType(outputDef.wdlType)
    }.toMap
  }
}
