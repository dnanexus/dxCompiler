package dx.translator.cwl

import dx.api.{DxApi, DxPath}
import dx.core.Constants
import dx.core.io.DxWorkerPaths
import dx.core.ir.RunSpec.{DefaultInstanceType, NoImage}
import dx.core.ir.{
  Application,
  BlockKind,
  Callable,
  CallableAttribute,
  DxName,
  ExecutableKindApplet,
  ExecutableKindWfCustomReorgOutputs,
  ExecutableKindWfFragment,
  ExecutableKindWfOutputs,
  InstanceTypeSelection,
  Level,
  Parameter,
  ParameterAttribute,
  SourceCode,
  Stage,
  StageInput,
  StageInputEmpty,
  StageInputStageLink,
  StageInputStatic,
  StageInputWorkflowLink,
  Type,
  Value,
  Workflow
}
import dx.core.languages.cwl.{
  CwlBlock,
  CwlBundle,
  CwlDxName,
  CwlUtils,
  DxHints,
  OptionalBlockInput,
  RequiredBlockInput,
  RequirementEvaluator,
  TargetParam
}
import dx.cwl.{
  CommandLineTool,
  CommandOutputParameter,
  CwlDirectory,
  CwlFile,
  CwlMulti,
  CwlOptional,
  CwlType,
  CwlValue,
  DirectoryValue,
  Evaluator,
  EvaluatorContext,
  ExpressionTool,
  FileValue,
  Hint,
  InputParameter,
  LoadListing,
  OutputParameter,
  PathValue,
  Process,
  Requirement,
  Runtime,
  StringValue,
  WorkflowInputParameter,
  WorkflowOutputParameter,
  WorkflowStep,
  WorkflowStepInput,
  Parameter => CwlParameter,
  Workflow => CwlWorkflow
}
import dx.translator.CallableAttributes.{DescriptionAttribute, TitleAttribute}
import dx.translator.ParameterAttributes.{HelpAttribute, LabelAttribute}
import dx.translator.{CustomReorgSettings, DxWorkflowAttrs, ReorgSettings, WorkflowTranslator}
import dx.util.protocols.{DxFileSource, DxFolderSource}
import dx.util.{FileSourceResolver, FileUtils, LocalFileSource, Logger}
import spray.json._

import java.nio.file.{Path, Paths}
import scala.util.Try

case class CwlSourceCode(source: Path, name: String = "") extends SourceCode {
  override val language: String = "cwl"
  override def toString: String = FileUtils.readFileContent(source)
  override def getDocContents: String = {
    def inner(jsv: JsValue, name: String): Option[JsValue] = {
      jsv match {
        case JsObject(fields) if fields.get("id").contains(JsString(s"${name}")) => Some(jsv)
        case JsObject(fields) if fields.contains("run") =>
          inner(fields.getOrElse("run", JsNull), name)
        case JsObject(fields) if fields.contains("steps") =>
          inner(fields.getOrElse("steps", JsArray.empty), name)
        case JsArray(items) =>
          val v = items.map(inner(_, name)).flatten
          v match {
            case Vector(i) => Some(i)
            case _         => None
          }
        case _ => None
      }
    }
    if (name.nonEmpty) {
      inner(toString.parseJson, name).getOrElse(JsString.empty).sortedPrint
    } else {
      toString.parseJson.sortedPrint
    }
  }
}

case class ProcessTranslator(cwlBundle: CwlBundle,
                             cwlSchemas: Option[JsValue],
                             locked: Boolean,
                             defaultRuntimeAttrs: Map[String, Value],
                             reorgAttrs: ReorgSettings,
                             perWorkflowAttrs: Map[String, DxWorkflowAttrs],
                             defaultScatterChunkSize: Int,
                             useManifests: Boolean,
                             instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
                             dxApi: DxApi = DxApi.get,
                             fileResolver: FileSourceResolver = FileSourceResolver.get,
                             logger: Logger = Logger.get) {

  private lazy val schemaStaticDependencies = cwlSchemas match {
    case Some(JsArray(schemas)) =>
      schemas.collect {
        case JsString(uri) if uri.startsWith(DxPath.DxUriPrefix) => uri
      }.toSet
    case None => Set.empty
    case Some(other) =>
      throw new Exception(s"invalid $$schemas value ${other}")
  }
  private lazy val cwlDefaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] =
    CwlUtils.fromIRValues(defaultRuntimeAttrs, isInput = true)

  private def translateParameterAttributes(
      param: CwlParameter,
      hintParameterAttrs: Map[String, Vector[ParameterAttribute]]
  ): Vector[ParameterAttribute] = {
    param.getName
      .map { name =>
        val metaAttrs =
          Vector(param.doc.map(HelpAttribute), param.label.map(LabelAttribute)).flatten
        val hintAttrs = hintParameterAttrs.getOrElse(name, Vector.empty)
        metaAttrs ++ hintAttrs
      }
      .getOrElse(Vector.empty)
  }

  private def translateCallableAttributes(
      process: Process,
      hintCallableAttributes: Map[String, Vector[CallableAttribute]]
  ): Vector[CallableAttribute] = {
    val metaAttrs =
      Vector(process.doc.map(DescriptionAttribute), process.label.map(TitleAttribute)).flatten
    val hintAttrs = hintCallableAttributes.getOrElse(process.name, Vector.empty)
    metaAttrs ++ hintAttrs
  }

  private def evaluateDefault(
      default: CwlValue,
      cwlType: CwlType,
      evaluator: Evaluator,
      evaluatorContext: EvaluatorContext = EvaluatorContext.empty
  ): Option[Value] = {
    // We need to deal with the case where a default value is a string and it's type is File or
    // Directory (or a multi-type that contains File/Directory). Only remote files can be specified
    // as default values. If the string represents a remote file, then it's treated as a file,
    // otherwise its coerced to a different type if possible (in the multi-type case), otherwise None.
    val result =
      try {
        evaluator.evaluate(default, cwlType, evaluatorContext) match {
          case (CwlMulti(types), s: StringValue)
              if types.map(CwlOptional.unwrapOptional).contains(CwlFile) =>
            if (s.value.startsWith(DxPath.DxUriPrefix)) {
              Some(s.coerceTo(CwlFile))
            } else if (types.size > 1) {
              Some(
                  s.coerceTo(CwlMulti(types.filter(t => CwlOptional.unwrapOptional(t) != CwlFile)))
              )
            } else {
              None
            }
          case (CwlMulti(types), s: StringValue)
              if types.map(CwlOptional.unwrapOptional).contains(CwlDirectory) =>
            if (s.value.startsWith(DxPath.DxUriPrefix)) {
              Some(s.coerceTo(CwlDirectory))
            } else if (types.size > 1) {
              Some(
                  s.coerceTo(
                      CwlMulti(types.filter(t => CwlOptional.unwrapOptional(t) != CwlDirectory))
                  )
              )
            } else {
              None
            }
          case (_, file: FileValue) if !CwlUtils.isDxPath(file)    => None
          case (_, dir: DirectoryValue) if !CwlUtils.isDxPath(dir) => None
          case (t, v)                                              => Some(v.coerceTo(t))
        }
      } catch {
        case _: Throwable => None
      }
    result.map {
      case (actualType, defaultValue) => CwlUtils.toIRValue(defaultValue, actualType)._2
    }
  }

  private case class CwlToolTranslator(
      tool: Process,
      isPrimary: Boolean,
      inheritedRequirements: Vector[Requirement] = Vector.empty,
      inheritedHints: Vector[Hint] = Vector.empty
  ) {
    private lazy val requirementEvaluator = RequirementEvaluator(
        inheritedRequirements ++ tool.requirements,
        inheritedHints ++ tool.hints,
        Map.empty,
        DxWorkerPaths.default,
        defaultRuntimeAttrs = cwlDefaultRuntimeAttrs,
        fileResolver = fileResolver
    )
    private lazy val dxHints = requirementEvaluator.getHintOfType[DxHints].map {
      case (dxHints: DxHints, _) => dxHints
    }
    private lazy val hintParameterAttrs: Map[String, Vector[ParameterAttribute]] = {
      dxHints
        .map(_.getParameterAttributes)
        .getOrElse(Map.empty)
    }
    private lazy val hintCallableAttrs: Map[String, Vector[CallableAttribute]] = {
      dxHints.map(_.getCallableAttributes).getOrElse(Map.empty)
    }

    def translateInput(input: InputParameter, evaluatorContext: EvaluatorContext): Parameter = {
      val irDefaultValue = input.default.flatMap(
          evaluateDefault(_, input.cwlType, requirementEvaluator.evaluator, evaluatorContext)
      )
      val dxName = CwlDxName.fromSourceName(input.name)
      val attrs = translateParameterAttributes(input, hintParameterAttrs)
      Parameter(dxName, CwlUtils.toIRType(input.cwlType), irDefaultValue, attrs)
    }

    def translateOutput(output: OutputParameter, evaluatorContext: EvaluatorContext): Parameter = {
      val irValue = output match {
        case cmdOutput: CommandOutputParameter =>
          cmdOutput.outputBinding
            .flatMap(_.outputEval.flatMap { cwlValue =>
              try {
                val (actualType, actualValue) =
                  requirementEvaluator.evaluator.evaluate(cwlValue,
                                                          output.cwlType,
                                                          evaluatorContext,
                                                          coerce = true)
                val (_, value) = CwlUtils.toIRValue(actualValue, actualType)
                Some(value)
              } catch {
                // the expression cannot be statically evaluated
                case _: Throwable => None
              }
            })
        case _ => None
      }
      val dxName = CwlDxName.fromSourceName(output.name)
      val attrs = translateParameterAttributes(output, hintParameterAttrs)
      Parameter(dxName, CwlUtils.toIRType(output.cwlType), irValue, attrs)
    }

    // TODO: If the source CWL has a DockerRequirement with a dockerLoad URI like
    //  dx://dxCompiler_playground:/glnexus_internal, rewrite it to dx://project-xxx:file-xxx to
    //  avoid a runtime lookup
    def apply: Application = {
      val name = tool.name
      // translate inputs and evaluate any static default values
      val inputParams = tool.inputs.collect {
        case i if i.id.isDefined => i.name -> i
      }.toMap
      val inputCtx = CwlUtils.createEvaluatorContext(runtime = Runtime.empty)
      // cwl applications always have a "target___" parameter
      val inputs = inputParams.values.map(i => translateInput(i, inputCtx)).toVector :+ TargetParam
      val defaults = inputParams.values.collect {
        case i if i.default.isDefined => i.name -> (i.cwlType, i.default.get)
      }.toMap
      // translate outputs and evaluate any static expressions
      val inputDir =
        tool.source.flatMap(src => Option(Paths.get(src).getParent)).getOrElse(Paths.get("."))
      val outputCtx = CwlUtils.createEvaluatorContext(env = defaults,
                                                      inputParameters = inputParams,
                                                      inputDir = inputDir,
                                                      fileResolver = fileResolver,
                                                      runtime = Runtime.empty)
      val outputs = tool.outputs.collect {
        case i if i.id.isDefined => translateOutput(i, outputCtx)
      }
      val docSource = tool.source.orElse(cwlBundle.primaryProcess.source) match {
        case Some(path) => Paths.get(path)
        case None       => throw new Exception(s"no source code for tool ${name}")
      }

      def extractPaths(p: PathValue): Vector[String] = {
        p match {
          case f: FileValue if f.location.isDefined || f.path.isDefined =>
            val fileUri = fileResolver.resolve(f.location.orElse(f.path).get) match {
              case DxFileSource(dxFile, _) => Vector(dxFile.asUri)
              case local: LocalFileSource =>
                logger.warning(
                    s"""InitialWorkDirRequirement in ${tool.name} or one of its ancestors 
                       |references local path ${local.canonicalPath}. This path must be
                       |already staged on the worker or in the Docker container, otherwise
                       |you must first upload the file and replace the location with a
                       |URI of the form 'dx://<project-id>:<file-id>'.""".stripMargin
                      .replaceAll("\n", " ")
                )
                Vector.empty
              case _ => Vector.empty
            }
            fileUri ++ f.secondaryFiles.flatMap(extractPaths)
          case d: DirectoryValue if d.location.isDefined =>
            fileResolver.resolveDirectory(d.location.get) match {
              case fs: DxFolderSource => Vector(fs.address)
              case local: LocalFileSource =>
                logger.warning(
                    s"""InitialWorkDirRequirement in ${tool.name} or one of its ancestors 
                       |references local path ${local.canonicalPath}. This path must be
                       |already staged on the worker or in the Docker container, otherwise
                       |you must first upload the directory and replace the location with a
                       |URI of the form 'dx://<project-id>:<folder>'.""".stripMargin
                      .replaceAll("\n", " ")
                )
                Vector.empty
              case _ => Vector.empty
            }
          case d: DirectoryValue => d.listing.flatMap(extractPaths)
          case _                 => Vector.empty
        }
      }
      val standAloneTask = CwlSourceCode(docSource, name)
      Application(
          name,
          inputs,
          outputs,
          // Try to determine the instance type from the process requirements -
          // if any of the expressions fail to evaluate statically, then the instance
          // type will be determined at runtime.
          requirementEvaluator.translateInstanceType(instanceTypeSelection),
          requirementEvaluator.translateContainer,
          ExecutableKindApplet,
          standAloneTask,
          translateCallableAttributes(tool, hintCallableAttrs),
          requirementEvaluator.translateApplicationRequirements,
          tags = Set("cwl"),
          // extract all paths that are hard-coded in the CWL
          staticFileDependencies =
            schemaStaticDependencies ++
              requirementEvaluator.translatePathDependencies.flatMap(extractPaths).toSet
      )
    }
  }

  private case class CwlWorkflowTranslator(
      wf: CwlWorkflow,
      availableDependencies: Map[String, Callable],
      workflowAttrs: Option[DxWorkflowAttrs],
      inheritedRequirements: Vector[Requirement] = Vector.empty,
      inheritedHints: Vector[Hint] = Vector.empty
  ) extends WorkflowTranslator(wf.name, availableDependencies, reorgAttrs, logger) {
    // Only the toplevel workflow may be unlocked, and only if it has
    // at least one step. This happens only if the user specifically
    // compiles it as "unlocked".
    protected lazy val isLocked: Boolean = {
      wf.steps.isEmpty || (cwlBundle.primaryProcess match {
        case wf2: CwlWorkflow => (wf.name != wf2.name) || locked
        case _                => true
      })
    }

    private lazy val evaluator: Evaluator = Evaluator.default

    private case class CwlCallEnv(env: Map[DxName, LinkedVar]) extends CallEnv(env) {
      override protected def create(env: Map[DxName, (Parameter, StageInput)]): CallEnv = {
        CwlCallEnv(env)
      }

      override def lookup(dxName: DxName): Option[(DxName, (Parameter, StageInput))] = {
        DxName.lookup(dxName, env)
      }
    }

    private object CwlCallEnv {
      def fromLinkedVars(lvars: Vector[LinkedVar]): CallEnv = {
        CwlCallEnv(lvars.map {
          case (parameter, stageInput) => parameter.name -> (parameter, stageInput)
        }.toMap)
      }
    }

    private lazy val docSource: Path =
      wf.source.orElse(cwlBundle.primaryProcess.source) match {
        case Some(path) => Paths.get(path)
        case None       => throw new Exception(s"no source code for tool ${wf.name}")
      }

    override protected def standAloneWorkflow: SourceCode = {
      CwlSourceCode(docSource, wf.name)
    }

    private def createWorkflowInput(input: WorkflowInputParameter): Parameter = {
      val dxName = CwlDxName.fromSourceName(input.name)
      val cwlType = input.cwlType
      val irType = CwlUtils.toIRType(cwlType)
      val default = input.default.flatMap(evaluateDefault(_, cwlType, evaluator))
      Parameter(dxName, irType, default, Vector.empty)
    }

    private def callInputToStageInput(callInput: Option[WorkflowStepInput],
                                      calleeParam: Parameter,
                                      env: CallEnv,
                                      locked: Boolean,
                                      callName: String): StageInput = {

      def lookup(key: DxName): StageInput = {
        // the key may be prefixed by a workflow or stage namespace - we start with the full name
        // and remove the namespaces successively (left to right) until we find a match
        env
          .lookup(key)
          .map {
            case (_, (_, stageInput)) => stageInput
          }
          .getOrElse(
              throw new Exception(
                  s"""|input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callName}>
                      |is missing from the environment. We don't have '${key}' in the environment.
                      |""".stripMargin.replaceAll("\n", " ")
              )
          )
      }

      def getDefaultInput(stepInput: WorkflowStepInput): StageInput = {
        Try {
          stepInput.default.map { d =>
            val cwlType = CwlUtils.fromIRType(calleeParam.dxType, isInput = true)
            val (_, cwlValue) = CwlUtils.toIRValue(d, cwlType)
            StageInputStatic(cwlValue)
          }
        }.toOption.flatten.getOrElse(StageInputEmpty)
      }

      callInput match {
        case None if Type.isOptional(calleeParam.dxType) => StageInputEmpty
        case None if calleeParam.defaultValue.isDefined  => StageInputEmpty
        case None if locked =>
          env.log()
          throw new Exception(
              s"""|input <${calleeParam.name}, ${calleeParam.dxType}> to call <${callName}>
                  |is unspecified. This is illegal in a locked workflow.""".stripMargin
                .replaceAll("\n", " ")
          )
        case None =>
          // the callee may not use the argument - defer until runtime
          StageInputEmpty
        case Some(inp) if inp.sources.size == 1 =>
          val dxName = CwlDxName.fromDecodedName(inp.sources.head.frag)
          try {
            lookup(dxName)
          } catch {
            case _: Throwable => getDefaultInput(inp)
          }
        case Some(inp) if inp.sources.isEmpty => getDefaultInput(inp)
        case Some(inp) =>
          throw new Exception(s"call input has multiple sources: ${inp}")
      }
    }

    private def getTarget(calleeName: String): String = {
      cwlBundle.processParents.get(calleeName) match {
        case None | Some(Vector()) => throw new Exception(s"process not found ${calleeName}")
        case Some(Vector(target))  => target
        case Some(target)          => s"${target.head}#${target.drop(1).mkString("/")}"
      }
    }

    private def translateCall(call: WorkflowStep, env: CallEnv, locked: Boolean): Stage = {
      val calleeName = call.run.name
      val callee: Callable = availableDependencies.getOrElse(
          calleeName,
          throw new Exception(
              s"""|Callable ${calleeName} should exist but is missing from the list of known 
                  |tasks/workflows ${availableDependencies.keys}|""".stripMargin
                .replaceAll("\n", " ")
          )
      )
      // Extract the input values/links from the environment
      val callInputParams: Map[DxName, WorkflowStepInput] = call.inputs.map { param =>
        CwlDxName.fromSourceName(param.name) -> param
      }.toMap
      val inputs: Vector[StageInput] = callee.inputVars.map {
        case param if param == TargetParam => StageInputStatic(Value.VString(getTarget(calleeName)))
        case param =>
          callInputToStageInput(callInputParams.get(param.name), param, env, locked, call.name)
      }
      Stage(call.name, getStage(), calleeName, inputs, callee.outputVars)
    }

    /**
      * Builds an applet to evaluate a WDL workflow fragment.
      *
      * @param wfName the workflow name
      * @param block the Block to translate into a WfFragment
      * @param blockPath keeps track of which block this fragment represents;
      *                   a top level block is a number. A sub-block of a top-level
      *                   block is a vector of two numbers, etc.
      * @param stepPath the workflow steps that led to this fragment.
      * @param env the environment
      * @return
      */
    private def translateWfFragment(wfName: String,
                                    block: CwlBlock,
                                    blockPath: Vector[Int],
                                    stepPath: Option[String],
                                    env: CallEnv): (Stage, Callable) = {
      val stageName = block.getName match {
        case None       => Constants.EvalStage
        case Some(name) => name
      }
      logger.trace(s"Compiling fragment <$stageName> as stage")
      logger
        .withTraceIfContainsKey("GenerateIR")
        .trace(
            s"""|block:
                |${block.prettyFormat}
                |""".stripMargin
        )

      // Get the applet inputs from the block inputs, which represent the closure required for the
      // block - all the variables defined earlier that are required to evaluate any expression.
      // Some referenced variables may be undefined because they are optional.
      val (inputParams, stageInputs) = block.inputs.flatMap {
        case RequiredBlockInput(name, _) =>
          env.lookup(name) match {
            case Some((fqn, (param, stageInput))) =>
              Some(param.copy(name = fqn), stageInput)
            case None =>
              throw new Exception(s"missing required input ${name}")
          }
        case OptionalBlockInput(name, _) =>
          env.lookup(name).map {
            case (fqn, (param, stageInput)) => (param.copy(name = fqn), stageInput)
          }
      }.unzip

      // The fragment runner can only handle a single call. If the block already has exactly one
      // call, then we translate it as an applet. If the block contains a scatter/conditional with
      // several calls, then we translate the inner block into a sub-workflow. Also determines the
      // name of the callable - we need to link with it when we get to the compile phase.
      val (calleeName, newStepPath) = block.kind match {
        case BlockKind.CallDirect =>
          throw new Exception(s"a direct call should not reach this stage")
        case BlockKind.CallFragment | BlockKind.ConditionalOneCall =>
          // a block with no nested sub-blocks, and a single call, or
          // a conditional with exactly one call in the sub-block
          (Some(block.target.run.name), None)
        case BlockKind.ScatterOneCall =>
          // a scatter with exactly one call in the sub-block
          val stepName = block.target.run.name
          val newScatterPath = stepPath.map(p => s"${p}.${stepName}").getOrElse(stepName)
          (Some(stepName), Some(newScatterPath))
        case _ =>
          throw new Exception(s"unexpected block ${block.prettyFormat}")
      }

      val scatterChunkSize: Option[Int] = newStepPath.map { sctPath =>
        val scatterChunkSize = workflowAttrs
          .flatMap { wfAttrs =>
            wfAttrs.scatters
              .flatMap {
                _.get(sctPath)
                  .orElse(wfAttrs.scatterDefaults)
                  .flatMap(scatterAttrs => scatterAttrs.chunkSize)
              }
          }
          .getOrElse(defaultScatterChunkSize)
        if (scatterChunkSize > Constants.JobsPerScatterLimit) {
          logger.warning(
              s"The number of jobs per scatter must be between 1-${Constants.JobsPerScatterLimit}"
          )
          Constants.JobsPerScatterLimit
        } else {
          scatterChunkSize
        }
      }

      val outputParams: Vector[Parameter] = block.outputs.map { param =>
        val irType = CwlUtils.toIRType(param.cwlType)
        Parameter(param.name, irType)
      }

      // create the type map that will be serialized in the applet's details
      val fqnDictTypes: Map[DxName, Type] = inputParams.map { param: Parameter =>
        param.name -> param.dxType
      }.toMap

      val standAloneFrag = CwlSourceCode(docSource, calleeName.get)
      val applet = Application(
          s"${wfName}_frag_${getStageId()}",
          inputParams,
          outputParams,
          DefaultInstanceType,
          NoImage,
          ExecutableKindWfFragment(calleeName, blockPath, fqnDictTypes, scatterChunkSize),
          standAloneFrag
      )

      (Stage(stageName, getStage(), applet.name, stageInputs, outputParams), applet)
    }

    private def createWorkflowStages(
        wfName: String,
        wfInputs: Vector[LinkedVar],
        blockPath: Vector[Int],
        subBlocks: Vector[CwlBlock],
        stepPath: Option[String],
        locked: Boolean
    ): (Vector[(Stage, Vector[Callable])], CallEnv) = {
      logger.trace(s"Assembling workflow backbone $wfName")

      val inputEnv: CallEnv = CwlCallEnv.fromLinkedVars(wfInputs)

      val logger2 = logger.withIncTraceIndent()
      logger2.trace(s"inputs: ${inputEnv.keys}")

      // link together all the stages into a linear workflow
      val (allStageInfo, stageEnv): (Vector[(Stage, Vector[Callable])], CallEnv) =
        subBlocks.zipWithIndex.foldLeft((Vector.empty[(Stage, Vector[Callable])], inputEnv)) {
          case ((stages, beforeEnv), (block: CwlBlock, blockNum: Int)) =>
            if (block.kind != BlockKind.CallDirect) {
              // A simple block that requires just one applet, OR
              // a complex block that needs a subworkflow
              val (stage, callable) =
                translateWfFragment(wfName, block, blockPath :+ blockNum, stepPath, beforeEnv)
              val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                case (env, param) =>
                  env.add(param.name, (param, StageInputStageLink(stage.dxStage, param)))
              }
              (stages :+ (stage, Vector(callable)), afterEnv)
            } else if (block.targetIsSimpleCall) {
              val step = block.target
              // The block contains exactly one call, with no extra variables. All the variables are
              // already in the environment, so there is no need to do any extra work. Compile
              // directly into a workflow stage.
              logger2.trace(s"Translating step ${step.name} as stage")
              val stage = translateCall(step, beforeEnv, locked)
              // Add bindings for the output variables. This allows later calls to refer
              // to these results.
              val afterEnv = stage.outputs.foldLeft(beforeEnv) {
                case (env, param: Parameter) =>
                  val fqn = param.name.pushDecodedNamespace(step.name)
                  val paramFqn = param.copy(name = fqn)
                  env.add(fqn, (paramFqn, StageInputStageLink(stage.dxStage, param)))
              }
              (stages :+ (stage, Vector.empty[Callable]), afterEnv)
            } else {
              throw new Exception(s"invalid DirectCall block ${block}")
            }
        }

      if (logger2.containsKey("GenerateIR")) {
        logger2.trace(s"stages for workflow $wfName = [")
        val logger3 = logger2.withTraceIfContainsKey("GenerateIR", indentInc = 1)
        allStageInfo.foreach {
          case (stage, _) =>
            logger3.trace(
                s"${stage.description}, ${stage.dxStage.id} -> callee=${stage.calleeName}"
            )
        }
        logger2.trace("]")
      }

      (allStageInfo, stageEnv)
    }

    /**
      * Build an applet + workflow stage for evaluating outputs. There are two reasons
      * to be build a special output section:
      * 1. Locked workflow: some of the workflow outputs are expressions.
      *    We need an extra applet+stage to evaluate them.
      * 2. Unlocked workflow: there are no workflow outputs, so we create
      *    them artificially with a separate stage that collects the outputs.
      * @param wfName the workflow name
      * @param outputs the outputs
      * @param env the environment
      * @return
      */
    private def createOutputStage(wfName: String,
                                  outputs: Vector[WorkflowOutputParameter],
                                  blockPath: Vector[Int],
                                  env: CallEnv): (Stage, Application) = {
      val paramNames: Set[DxName] = outputs.flatMap { param =>
        if (param.sources.isEmpty) {
          throw new Exception(
              s"there must be at least one 'outputSource' for parameter ${param.name}"
          )
        }
        param.sources.map(source => CwlDxName.fromDecodedName(source.frag))
      }.toSet
      logger.trace(s"paramNames: ${paramNames}")
      val (applicationInputs, stageInputs) = paramNames.map { name =>
        env.lookup(name) match {
          case Some((_, (param, stageInput))) => (param, stageInput)
          case None =>
            throw new Exception(s"parameter ${name} missing from CallEnv")
        }
      }.unzip

      // build definitions of the output variables - if the expression can be evaluated,
      // set the values as the parameter's default
      val outputParams: Vector[Parameter] = outputs.map { param =>
        val dxName = CwlDxName.fromSourceName(param.name)
        val irType = CwlUtils.toIRType(param.cwlType)
        Parameter(dxName, irType, None)
      }

      // Determine kind of application. If a custom reorg app is used and this is a top-level
      // workflow (custom reorg applet doesn't apply to locked workflows), add an output
      // variable for reorg status.
      val (applicationKind, updatedOutputVars) = reorgAttrs match {
        case CustomReorgSettings(_, _, true) if !isLocked =>
          val updatedOutputVars = outputParams :+ Parameter(
              Constants.ReorgStatus,
              Type.TString,
              Some(Value.VString(Constants.ReorgStatusCompleted))
          )
          (ExecutableKindWfCustomReorgOutputs, updatedOutputVars)
        case _ =>
          (ExecutableKindWfOutputs(blockPath), outputParams)
      }
      val application = Application(
          s"${wfName}_${Constants.OutputStage}",
          applicationInputs.toVector,
          updatedOutputVars,
          DefaultInstanceType,
          NoImage,
          applicationKind,
          standAloneWorkflow
      )
      val stage = Stage(
          Constants.OutputStage,
          getStage(Some(Constants.OutputStage)),
          application.name,
          stageInputs.toVector,
          updatedOutputVars
      )
      (stage, application)
    }

    /**
      * Compile a locked workflow. This is called at the top level for locked workflows,
      * and it is always called for nested workflows regardless of whether the top level
      * is locked.
      * @param wfName workflow name
      * @param inputs formal workflow inputs
      * @param outputs workflow outputs
      * @param blockPath the path to the current (sub)workflow, as a vector of block indices
      * @param subBlocks the sub-blocks of the current block
      * @param stepPath the workflow steps that led to this fragment.
      * @param level the workflow level
      * @return
      */
    def translateWorkflowLocked(
        wfName: String,
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        blockPath: Vector[Int],
        subBlocks: Vector[CwlBlock],
        stepPath: Option[String],
        level: Level.Value
    ): (Workflow, Vector[Callable], Vector[(Parameter, StageInput)]) = {
      val wfInputParams = inputs.map(createWorkflowInput)
      val wfInputLinks: Vector[LinkedVar] = wfInputParams.map(p => (p, StageInputWorkflowLink(p)))
      val loadAnyInputs =
        inputs.iterator.exists(i => i.loadContents || i.loadListing != LoadListing.No)
      val (backboneInputs, commonStageInfo) = if (useManifests || loadAnyInputs) {
        // If we are using manifests, we need an initial applet to merge multiple
        // manifests into a single manifest.
        val commonStageInputs = wfInputParams.map(p => StageInputWorkflowLink(p))
        val inputOutputs: Vector[Parameter] = inputs.map { i =>
          val dxName = CwlDxName.fromSourceName(i.name)
          Parameter(dxName, CwlUtils.toIRType(i.cwlType))
        }
        val (commonStage, commonApplet) =
          createCommonApplet(wfName, wfInputParams, commonStageInputs, inputOutputs, blockPath)
        val fauxWfInputs: Vector[LinkedVar] = commonStage.outputs.map { param =>
          val link = StageInputStageLink(commonStage.dxStage, param)
          (param, link)
        }
        (fauxWfInputs, Vector((commonStage, Vector(commonApplet))))
      } else {
        (wfInputLinks, Vector.empty)
      }

      // translate the Block(s) into workflow stages
      val (backboneStageInfo, env) = createWorkflowStages(
          wfName,
          backboneInputs,
          blockPath,
          subBlocks,
          stepPath,
          locked = true
      )
      val (stages, auxCallables) = (commonStageInfo ++ backboneStageInfo).unzip

      // We need a common output stage for any of these reasons:
      // 1. we need to build an output manifest
      // 2. there are no workflow steps (unlikely, but there is at least one conformance
      // test where this is the case), then the workflow will consist of only the
      // output applet.
      val useOutputStage = useManifests || backboneStageInfo.isEmpty || {
        // 3. there are outputs that need to be merged
        // 4. we are "downcasting" an output, e.g. the applet output is `Any` but
        // the workflow output is `File`
        // 5. an output is used directly as an input
        // For example, in the small workflow below, 'lane' is used in such a manner.
        //
        // workflow inner {
        //   input {
        //      String lane
        //   }
        //   output {
        //      String blah = lane
        //   }
        // }
        //
        // In locked workflows, it is illegal to access a workflow input directly from
        // a workflow output. It is only allowed to access a stage input/output.
        val inputNames = inputs.map(_.name).toSet
        outputs.exists { out =>
          out.sources.size > 1 ||
          out.linkMerge.isDefined ||
          out.pickValue.isDefined ||
          out.sources.exists(i => inputNames.contains(i.frag)) ||
          out.sources.headOption.exists { sourceId =>
            env.get(CwlDxName.fromDecodedName(sourceId.frag)).map(_._1).exists { stageParam =>
              try {
                CwlUtils.requiresDowncast(CwlUtils.fromIRType(stageParam.dxType, isInput = false),
                                          out.cwlType)
              } catch {
                case _: Throwable =>
                  // There is an edge case in CWL where a stage output may be optional due to the
                  // stage having a `when` clause, but it is referenced by a non-optional workflow
                  // output. The expectation is that there will be a runtime error when the condition
                  // is false (see https://github.com/common-workflow-language/cwl-v1.2/issues/146).
                  // We handle that in the output applet.
                  true
              }
            }
          }
        }
      }

      val (wfOutputs, finalStages, finalCallables) = if (useOutputStage) {
        val (outputStage, outputApplet) = createOutputStage(wfName, outputs, blockPath, env)
        val wfOutputs = outputStage.outputs.map { param =>
          (param, StageInputStageLink(outputStage.dxStage, param))
        }
        (wfOutputs, stages :+ outputStage, auxCallables.flatten :+ outputApplet)
      } else {
        val wfOutputs = outputs.map { out =>
          val dxName = CwlDxName.fromSourceName(out.name)
          val outputStage = env.get(dxName).map(_._2).getOrElse {
            // we know there is only one output source otherwise we'd be
            // using an output stage
            val outputSource = CwlDxName.fromDecodedName(out.sources.head.frag)
            env
              .lookup(outputSource)
              .map(_._2._2)
              .getOrElse(
                  throw new Exception(
                      s"output source ${outputSource} missing from environment ${env}"
                  )
              )
          }
          val irType = CwlUtils.toIRType(out.cwlType)
          (Parameter(dxName, irType), outputStage)
        }
        (wfOutputs, stages, auxCallables.flatten)
      }

      (Workflow(wfName,
                wfInputLinks,
                wfOutputs,
                finalStages,
                standAloneWorkflow,
                locked = true,
                level),
       finalCallables,
       wfOutputs)
    }

    private def translateTopWorkflowLocked(
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        subBlocks: Vector[CwlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      translateWorkflowLocked(wf.name, inputs, outputs, Vector.empty, subBlocks, None, Level.Top)
    }

    private def translateTopWorkflowUnlocked(
        inputs: Vector[WorkflowInputParameter],
        outputs: Vector[WorkflowOutputParameter],
        subBlocks: Vector[CwlBlock]
    ): (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      // Create a special applet+stage for the inputs. This is a substitute for workflow inputs. We
      // now call the workflow inputs, "fauxWfInputs" since they are references to the outputs of
      // this first applet.
      val commonAppletInputs: Vector[Parameter] =
        inputs.map(input => createWorkflowInput(input))
      val commonStageInputs: Vector[StageInput] = inputs.map(_ => StageInputEmpty)
      val (commonStg, commonApplet) =
        createCommonApplet(wf.name, commonAppletInputs, commonStageInputs, commonAppletInputs)
      val fauxWfInputs: Vector[LinkedVar] = commonStg.outputs.map { param =>
        val stageInput = StageInputStageLink(commonStg.dxStage, param)
        (param, stageInput)
      }

      val (allStageInfo, env) =
        createWorkflowStages(wf.name, fauxWfInputs, Vector.empty, subBlocks, None, locked = false)
      val (stages, auxCallables) = allStageInfo.unzip

      // convert the outputs into an applet+stage
      val (outputStage, outputApplet) =
        createOutputStage(wf.name, outputs, Vector.empty, env)

      val wfInputs = commonAppletInputs.map(param => (param, StageInputEmpty))
      val wfOutputs =
        outputStage.outputs.map(param => (param, StageInputStageLink(outputStage.dxStage, param)))

      val irwf = Workflow(
          wf.name,
          wfInputs,
          wfOutputs,
          commonStg +: stages :+ outputStage,
          standAloneWorkflow,
          locked = false,
          Level.Top
      )
      (irwf, commonApplet +: auxCallables.flatten :+ outputApplet, wfOutputs)
    }

    def translate: (Workflow, Vector[Callable], Vector[LinkedVar]) = {
      logger.trace(s"Translating workflow ${wf.name}")
      // Create a stage per workflow step - steps with conditional and/or scatter
      // require helper applets.
      val subBlocks = CwlBlock.createBlocks(wf)
      if (isLocked) {
        translateTopWorkflowLocked(wf.inputs, wf.outputs, subBlocks)
      } else {
        translateTopWorkflowUnlocked(wf.inputs, wf.outputs, subBlocks)
      }
    }
  }

  def translateProcess(process: Process,
                       availableDependencies: Map[String, Callable],
                       isPrimary: Boolean): Vector[Callable] = {
    process match {
      case tool: CommandLineTool =>
        val toolTranslator =
          CwlToolTranslator(tool,
                            isPrimary = isPrimary,
                            cwlBundle.requirements.getOrElse(tool.name, Vector.empty),
                            cwlBundle.hints.getOrElse(tool.name, Vector.empty))
        Vector(toolTranslator.apply)
      case tool: ExpressionTool =>
        val toolTranslator = CwlToolTranslator(
            tool,
            isPrimary = isPrimary,
            cwlBundle.requirements.getOrElse(tool.name, Vector.empty),
            cwlBundle.hints.getOrElse(tool.name, Vector.empty)
        )
        Vector(toolTranslator.apply)
      case wf: CwlWorkflow =>
        val wfAttrs = perWorkflowAttrs.get(wf.name)
        val wfTranslator = CwlWorkflowTranslator(
            wf,
            availableDependencies,
            wfAttrs,
            cwlBundle.requirements.getOrElse(wf.name, Vector.empty),
            cwlBundle.hints.getOrElse(wf.name, Vector.empty)
        )
        wfTranslator.apply
      case _ =>
        throw new Exception(s"Process type ${process} not supported")
    }
  }
}
