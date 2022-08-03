package dx.executor.cwl

import dx.api.{DxExecution, DxWorkflow}
import dx.core.Constants
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.core.ir.{Block, DxName, ExecutableLink, ParameterLink, Type, TypeSerde, Value, ValueSerde}
import dx.core.languages.Language
import dx.core.languages.cwl.{
  CwlBlock,
  CwlDxName,
  CwlUtils,
  DxHintSchema,
  OptionalBlockInput,
  RequiredBlockInput,
  RequirementEvaluator,
  Target,
  TargetParam
}
import dx.cwl.{
  ArrayValue,
  BooleanValue,
  CwlAny,
  CwlArray,
  CwlBoolean,
  CwlMulti,
  CwlOptional,
  CwlRecord,
  CwlType,
  CwlValue,
  CwlValueContext,
  DirectoryValue,
  Evaluator,
  EvaluatorContext,
  FileValue,
  HintUtils,
  Identifiable,
  Identifier,
  InputParameter,
  LinkMergeMethod,
  LoadListing,
  Loadable,
  NullValue,
  ObjectValue,
  Parameter,
  Parser,
  PickValueMethod,
  Process,
  ScatterMethod,
  Workflow,
  WorkflowInputParameter,
  WorkflowOutputParameter,
  WorkflowStepInput
}
import dx.executor.{JobMeta, WorkflowExecutor}
import dx.util.{AddressableFileSource, FileSource, LocalFileSource, TraceLevel}
import spray.json._

import java.net.URI
import scala.annotation.tailrec

object CwlWorkflowExecutor {
  private val nameRegex = "(.+)_frag_stage-.+".r

  def create(jobMeta: JobMeta, separateOutputs: Boolean): CwlWorkflowExecutor = {
    // when parsing a packed workflow as a String, we need to use a baseuri -
    // it doesn't matter what it is
    val parser = Parser.create(Some(URI.create("file:/null")), hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(jobMeta.sourceCode) match {
      case (version, Some("Workflow")) if Language.parse(version) == Language.CwlV1_2 =>
        ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL Workflow document of a supported version
               |${jobMeta.sourceCode}""".stripMargin
        )
    }
    val wfName = jobMeta.getExecutableDetail(Constants.OriginalName) match {
      case Some(JsString(nameRegex(name))) => name
      case Some(JsString(name))            => name
      case _                               => throw new Exception("missing executable name")
    }
    val parserResult = parser.parseString(jobMeta.sourceCode)
    // A CWL workflow may contain nested workflows, and we may be executing the top-level workflow
    // or one nested at any level. Here we find all (sub-)workflows with a matching name and
    // simplify them, which makes them comparable. This is necessary because the same workflow may
    // be imported multiple times, and as long as two workflows are identical we don't care which
    // version we execute. We also track the full id of the workflow because we will need this to
    // fully specify the target step when launching the target applet.
    val (workflow, parentId) = parserResult.document.values
      .foldLeft(Map.empty[Process, Identifier]) {
        case (accu, wf: Workflow) if wf.name == wfName =>
          accu + (CwlUtils.simplifyProcess(wf) -> wf.id.get)
        case (accu, _) => accu
      }
      .toVector match {
      case Vector((wf: Workflow, id: Identifier)) => (wf, id)
      case Vector() =>
        parserResult.mainProcess match {
          case Some(wf: Workflow) => (wf, wf.id.get)
          case _ =>
            throw new Exception(s"expected CWL document to contain a workflow named ${wfName}")
        }
      case _ =>
        throw new Exception(
            s"CWL document contains multiple workflows with name ${wfName}"
        )
    }
    CwlWorkflowExecutor(workflow, parentId, jobMeta, separateOutputs)
  }
}

case class CwlWorkflowExecutor(workflow: Workflow,
                               parentId: Identifier,
                               jobMeta: JobMeta,
                               separateOutputs: Boolean)
    extends WorkflowExecutor[CwlBlock](jobMeta, separateOutputs) {
  private val logger = jobMeta.logger
  private lazy val parent = CwlDxName.fromDecodedName(parentId.frag).getDecodedParts

  override val executorName: String = "dxExecutorCwl"

  override protected lazy val typeAliases: Map[String, Type.TSchema] = {
    HintUtils.getSchemaDefs(workflow.requirements).collect {
      case (name, schema: CwlRecord) => name -> CwlUtils.toIRSchema(schema)
    }
  }

  // gets the listing for a folder - if the folder does not exist, returns an empty listing
  private def getFolderListing(uri: String, recursive: Boolean): Vector[FileSource] = {
    val fs =
      try {
        jobMeta.fileResolver.resolveDirectory(uri)
      } catch {
        case _: Throwable => return Vector()
      }
    if (!fs.isDirectory) {
      throw new Exception(s"Directory-type input is not a directory: ${fs}")
    }
    if (fs.exists && fs.isListable) {
      fs.listing(recursive)
    } else {
      Vector()
    }
  }

  // creates a listing of VFile/VFolders from a nested listing of FileSources
  private def createListing(value: Vector[FileSource]): Vector[PathValue] = {
    value.map {
      case fs: AddressableFileSource if fs.isDirectory && fs.isListable =>
        val items = createListing(fs.listing())
        VFolder(fs.address, listing = Some(items))
      case fs: AddressableFileSource if fs.isDirectory => VFolder(fs.address)
      case fs: LocalFileSource                         => VFile(fs.address)
      case fs: AddressableFileSource                   => VFile(fs.address)
    }
  }

  override protected def evaluateInputs(
      jobInputs: Map[DxName, (Type, Value)]
  ): Map[DxName, (Type, Value)] = {
    // this might be the input for the entire workflow or just a subblock
    val inputs = jobMeta.blockPath match {
      case Vector() =>
        val workflowInputs: Map[DxName, WorkflowInputParameter] =
          workflow.inputs.map(param => CwlDxName.fromSourceName(param.name) -> param).toMap
        if (logger.isVerbose) {
          logger.trace(
              s"""input parameters:
                 |${workflowInputs
                   .map {
                     case (name, inp) =>
                       s"  ${CwlUtils.prettyFormatType(inp.cwlType)} ${name}"
                   }
                   .mkString("\n")}""".stripMargin
          )
        }
        jobInputs.map {
          case (name, (irType, file: VFile))
              if workflowInputs(name).loadContents && file.contents.isEmpty =>
            name -> (irType, file
              .copy(contents = Some(jobMeta.fileResolver.resolve(file.uri).readString)))
          case (name, (irType, dir: VFolder))
              if workflowInputs(name).loadListing != LoadListing.No && dir.listing.isEmpty =>
            val listing = createListing(
                getFolderListing(dir.uri,
                                 recursive = workflowInputs(name).loadListing == LoadListing.Deep)
            )
            name -> (irType, dir.copy(listing = Some(listing)))
          case i => i
        }
      case path =>
        val block: CwlBlock = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), path)
        val inputTypes = block.inputs.map {
          case RequiredBlockInput(name, cwlType, param: Parameter with Loadable) =>
            name -> (CwlUtils.toIRType(cwlType), param.loadContents, param.loadListing)
          case RequiredBlockInput(name, cwlType, param) =>
            name -> (CwlUtils.toIRType(cwlType), false, LoadListing.No)
          case OptionalBlockInput(name, cwlType, param: Parameter with Loadable) =>
            val irType = Type.ensureOptional(CwlUtils.toIRType(cwlType))
            name -> (irType, param.loadContents, param.loadListing)
          case OptionalBlockInput(name, cwlType, param) =>
            name -> (Type.ensureOptional(CwlUtils.toIRType(cwlType)), false, LoadListing.No)
        }.toMap
        if (logger.isVerbose) {
          logger.trace(
              s"""input parameters:
                 |${inputTypes
                   .map {
                     case (name, (irType, true, _)) =>
                       s"  ${TypeSerde.toString(irType)} ${name}"
                     case (name, (irType, _, loadListing)) if loadListing != LoadListing.No =>
                       s"  ${TypeSerde.toString(irType)} ${name} (loadListing=${loadListing})"
                     case (name, (irType, _, _)) =>
                       s"  ${TypeSerde.toString(irType)} ${name} (loadContents=true)"
                   }
                   .mkString("\n")}""".stripMargin
          )
        }
        jobInputs.collect {
          case (name, (_, v)) if inputTypes.contains(name) =>
            // coerce the input value to the target type
            val (irType, loadContents, loadListing) = inputTypes(name)
            val irValue = coerceTo(v, irType) match {
              case file: VFile if loadContents && file.contents.isEmpty =>
                file.copy(contents = Some(jobMeta.fileResolver.resolve(file.uri).readString))
              case dir: VFolder if loadListing != LoadListing.No && dir.listing.isEmpty =>
                val listing = createListing(
                    getFolderListing(dir.uri, recursive = loadListing == LoadListing.Deep)
                )
                dir.copy(listing = Some(listing))
              case value => value
            }
            name -> (irType, irValue)
          case i => i
        }
    }
    if (logger.isVerbose) {
      logger.trace(
          s"""input values:
             |${inputs
               .map {
                 case (name, (t, v)) =>
                   s"${TypeSerde.toString(t)} ${name} = ${ValueSerde.toString(v, verbose = true)}}"
               }
               .mkString("\n")}""".stripMargin
      )
    }
    inputs
  }

  case class CwlEnv[T](env: Map[DxName, T]) {
    def lookup(name: DxName): Option[T] = {
      env.get(name).orElse(env.get(name.dropNamespaces()))
    }

    def apply(name: DxName): T = {
      lookup(name).get
    }

    def update(withEnv: Map[DxName, T]): CwlEnv[T] = {
      CwlEnv(env ++ withEnv)
    }
  }

  override protected def evaluateOutputs(jobInputs: Map[DxName, (Type, Value)],
                                         addReorgStatus: Boolean): Map[DxName, (Type, Value)] = {
    // This might be the output for the entire workflow or just a subblock.
    val outputParams = jobMeta.blockPath match {
      case Vector() =>
        workflow.outputs.map { param =>
          (CwlDxName.fromSourceName(param.name), param, param.cwlType)
        }
      case path =>
        Block.getSubBlockAt(CwlBlock.createBlocks(workflow), path).outputs.map { out =>
          (out.name, out.source, out.cwlType)
        }
    }
    if (logger.isVerbose) {
      logger.trace(s"Evaluating output parameters:\n  ${outputParams.mkString("\n  ")}")
    }
    val env = CwlEnv(jobInputs)
    val irOutputs: Map[DxName, (Type, Value)] = outputParams.map {
      case (dxName, param: WorkflowOutputParameter, cwlType) if param.sources.nonEmpty =>
        val sourceValues = param.sources.map(src =>
          env.lookup(CwlDxName.fromDecodedName(src.frag)).getOrElse((TMulti.Any, VNull))
        )
        val irValue = if (param.linkMerge.nonEmpty || param.pickValue.nonEmpty) {
          val mergedValues = (sourceValues, param.linkMerge) match {
            case (_, Some(LinkMergeMethod.MergeFlattened)) =>
              sourceValues.flatMap {
                case (TArray(itemType, _), VArray(items)) => items.map(i => (itemType, i))
                case (t, value)                           => Vector((t, value))
              }
            case (Vector((t, v)), Some(LinkMergeMethod.MergeNested)) =>
              Vector((TArray(t), VArray(v)))
            case (Vector((TArray(itemType, _), VArray(items))), None) =>
              items.map(i => (itemType, i))
            case (Vector(other), None) =>
              throw new Exception(
                  s"""when there is a single output source and no linkMerge method, the output 
                     |source must be an array, not ${other}""".stripMargin.replaceAll("\n", " ")
              )
            case _ => sourceValues
          }
          if (param.pickValue.nonEmpty) {
            val nonNull = mergedValues.map(_._2).filterNot(_ == VNull)
            param.pickValue.get match {
              case PickValueMethod.FirstNonNull =>
                nonNull.headOption.getOrElse(
                    throw new Exception(s"all source values are null for parameter ${param.name}")
                )
              case PickValueMethod.TheOnlyNonNull if nonNull.size == 1 => nonNull.head
              case PickValueMethod.TheOnlyNonNull =>
                throw new Exception(
                    s"there is not exactly one non-null value for parameter ${param.name}"
                )
              case PickValueMethod.AllNonNull => VArray(nonNull)
            }
          } else {
            VArray(mergedValues.map(_._2))
          }
        } else if (sourceValues.size == 1) {
          sourceValues.head._2
        } else {
          VArray(sourceValues.map(_._2))
        }
        val irType = CwlUtils.toIRType(cwlType)
        dxName -> (irType, coerceTo(irValue, irType))
      case (dxName, _, cwlType) if jobInputs.contains(dxName) =>
        val irType = CwlUtils.toIRType(cwlType)
        dxName -> (irType, coerceTo(jobInputs(dxName)._2, irType))
      case (dxName, _, cwlType) if CwlOptional.isOptional(cwlType) =>
        dxName -> (CwlUtils.toIRType(cwlType), VNull)
      case (dxName, _, _) =>
        throw new Exception(s"missing required output ${dxName}")
    }.toMap
    if (addReorgStatus) {
      irOutputs + (Constants.ReorgStatus -> (TString, VString(Constants.ReorgStatusCompleted)))
    } else {
      irOutputs
    }
  }

  case class CwlBlockContext(block: CwlBlock, cwlEnv: CwlEnv[(CwlType, CwlValue)])
      extends BlockContext {
    private val step = block.target
    private lazy val runInputs: Map[DxName, InputParameter] = step.run.inputs.map { i =>
      CwlDxName.fromSourceName(i.name) -> i
    }.toMap
    private lazy val scatterNames = step.scatter.map(src => CwlDxName.fromSourceName(src.name))
    private lazy val eval = Evaluator.create(block.targetRequirements, block.targetHints)

    override lazy val env: Map[DxName, (Type, Value)] = CwlUtils.toIR(cwlEnv.env)

    // Each step input has zero or more sources. If there is no source, then the default value is
    // used or the value is null. Otherwise, linkMerge and pickValue are applied if specified.
    private def evaluateStepInput(
        stepInput: WorkflowStepInput,
        fullEnv: CwlEnv[(CwlType, CwlValue)]
    ): (DxName, (WorkflowStepInput, (CwlType, CwlValue))) = {
      val dxName = CwlDxName.fromSourceName(stepInput.name)
      val (cwlType, cwlValue) = Option
        .when(stepInput.sources.nonEmpty)({
          def itemsToArray(items: Vector[(CwlType, CwlValue)]): (CwlType, CwlValue) = {
            val (types, values) = items.unzip
            (CwlArray(CwlType.flatten(types.distinct)), ArrayValue(values))
          }
          val sourceValues =
            stepInput.sources.map(src => fullEnv(CwlDxName.fromDecodedName(src.frag)))
          val (cwlType, cwlValue) =
            if (stepInput.linkMerge.nonEmpty || stepInput.pickValue.nonEmpty) {
              val mergedValues = (sourceValues, stepInput.linkMerge) match {
                case (_, Some(LinkMergeMethod.MergeFlattened)) =>
                  sourceValues.map(s => (CwlOptional.unwrapOptional(s._1), s._2)).flatMap {
                    case (array: CwlArray, ArrayValue(items)) => items.map(i => (array.itemType, i))
                    case (t, value)                           => Vector((t, value))
                  }
                case (Vector((t, v)), Some(LinkMergeMethod.MergeNested)) =>
                  Vector((CwlArray(t), ArrayValue(Vector(v))))
                case (Vector((array: CwlArray, ArrayValue(items))), None) =>
                  items.map(i => (array.itemType, i))
                case (Vector(other), None) =>
                  throw new Exception(
                      s"""when there is a single output source and no linkMerge method, the output 
                         |source must be an array, not ${other}""".stripMargin
                        .replaceAll("\n", " ")
                  )
                case _ => sourceValues
              }
              if (stepInput.pickValue.nonEmpty) {
                val nonNull = mergedValues.filterNot(_._2 == NullValue)
                stepInput.pickValue.get match {
                  case PickValueMethod.FirstNonNull =>
                    nonNull.headOption.getOrElse(
                        throw new Exception(
                            s"all source values are null for parameter ${stepInput.name}"
                        )
                    )
                  case PickValueMethod.TheOnlyNonNull if nonNull.size == 1 => nonNull.head
                  case PickValueMethod.TheOnlyNonNull =>
                    throw new Exception(
                        s"there is not exactly one non-null value for parameter ${stepInput.name}"
                    )
                  case PickValueMethod.AllNonNull => itemsToArray(nonNull)
                }
              } else if (mergedValues.size == 1) {
                mergedValues.head
              } else {
                itemsToArray(mergedValues)
              }
            } else if (sourceValues.size == 1) {
              sourceValues.head
            } else {
              itemsToArray(sourceValues)
            }
          (cwlType, if (cwlValue == NullValue && stepInput.default.isDefined) {
            stepInput.default.get
          } else {
            cwlValue
          })
        })
        .map {
          case (cwlType, NullValue) if stepInput.default.isDefined =>
            // use the default if the value produced by source is null
            stepInput.default.get.coerceTo(cwlType)
          case other => other
        }
        .orElse(Option.when(stepInput.default.isDefined) {
          // use the default if there are no sources - infer the type from the value
          val cwlValue = stepInput.default.get
          val cwlType = CwlValue.inferType(cwlValue, CwlValueContext.Input)
          (cwlType, cwlValue)
        })
        .getOrElse((CwlOptional(CwlAny), NullValue))
      dxName -> (stepInput, (cwlType, cwlValue))
    }
    private def createEvalInputs(
        inputs: Map[WorkflowStepInput, (CwlType, CwlValue)]
    ): ObjectValue = {
      EvaluatorContext.createInputs(inputs.map {
        case (param, tv) =>
          val p: Identifiable with Loadable = param
          p -> tv
      }, fileResolver = jobMeta.fileResolver)
    }

    private def evaluateCallStepInputs(
        inputs: Option[Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))]] = None,
        isConditional: Boolean = false
    ): Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))] = {
      // Evaluate all the step inputs for a call. There may be step inputs that are not passed to
      // the callee but are referred to in expressions of other step inputs' valueFrom fields. Note
      // that the step input's type here might not match the callee's input type if valueFrom is
      // specified.
      val stepInputs = inputs.getOrElse(step.inputs.map(evaluateStepInput(_, cwlEnv)).toMap)
      // evaluate valueFrom expressions - create the evaluator lazily since it might not be needed
      lazy val evalInputs: ObjectValue = createEvalInputs(stepInputs.values.toMap)
      val stepInputValues = stepInputs.map {
        case (dxName, (stepInput, (sourceType, sourceValue))) =>
          if (stepInput.valueFrom.isDefined) {
            try {
              val cwlType = runInputs.get(dxName).map(_.cwlType).getOrElse(sourceType)
              val selfValue = EvaluatorContext.finalizeInputValue(
                  sourceValue,
                  sourceType,
                  stepInput,
                  fileResolver = jobMeta.fileResolver
              )
              dxName -> (stepInput, eval.evaluate(stepInput.valueFrom.get,
                                                  cwlType,
                                                  EvaluatorContext(selfValue, evalInputs),
                                                  coerce = true))
            } catch {
              case cause: Throwable =>
                throw new Exception(
                    s"error evaluating stepInput ${stepInput.name} valueFrom ${stepInput.valueFrom.get}",
                    cause
                )
            }
          } else {
            val cwlType = runInputs
              .get(dxName)
              .map {
                case inp if isConditional && CwlOptional.isOptional(sourceType) =>
                  // if the source type is optional, it may be used in the 'when' expression of a
                  // conditional, so we force the dest type to be optional
                  CwlOptional.ensureOptional(inp.cwlType)
                case inp => inp.cwlType
              }
              .getOrElse(sourceType)
            val (callType, callValue) =
              try {
                cwlType match {
                  case CwlAny =>
                    // if the process input has type Any but the step input is a more specific type,
                    // the call input type needs to be Any so that the value will be serialized as a
                    // hash
                    (CwlAny, sourceValue)
                  case _: CwlMulti =>
                    // if the process input has multiple allowed types but the step input is a more
                    // specific type, coerce the value but leave the call input type as CwlMulti
                    // so that the value will be serialized as a hash
                    val (_, coercedValue) = sourceValue.coerceTo(cwlType)
                    (cwlType, coercedValue)
                  case _ => sourceValue.coerceTo(cwlType)
                }
              } catch {
                case cause: Throwable =>
                  throw new Exception(
                      s"""effective type ${sourceType} of stepInput ${stepInput.name} value 
                         |${sourceValue} to parameter ${dxName} is not coercible to expected type 
                         |${cwlType}""".stripMargin.replaceAll("\n", " "),
                      cause
                  )
              }
            dxName -> (stepInput, (callType, callValue))
          }
      }
      if (logger.isVerbose) {
        logger.trace(s"Workflow step ${step.name} calling process ${step.run.name}")
        logger.traceLimited(s"Inputs:\n  ${stepInputValues.mkString("\n  ")}")
      }
      stepInputValues
    }

    /*
     * Format an identifier as the target workflow step to execute.
     * The block target id is in the format current_wf/step, and the parent is
     * the path from top-level wf to the current wf as a vector of string
     * Here the target string (which will be passed to task executor) will be formatted
     * as 'top_wf#path/to/current_wf/step'
     * For example, if the top-level workflow is 'top_wf' and we want
     * to run a step in a sub-workflow, it might be 'top_wf#step1/nested/step2'.
     */
    private lazy val target: String = {
      val targetId = block.target.id.get
      val dxName = CwlDxName.fromDecodedName(targetId.frag)
      val fullDxName = dxName.pushDecodedNamespaces(parent.dropRight(1))
      if (fullDxName.numParts == 1) {
        targetId.name
      } else {
        val (top_wf, step) = fullDxName.popDecodedNamespace()
        s"${top_wf}#${step.toString}"
      }
    }

    private def launchCall(
        callInputs: Map[DxName, (CwlType, CwlValue)],
        nameDetail: Option[String],
        folder: Option[String]
    ): (DxExecution, ExecutableLink, String) = {
      logger.traceLimited(
          s"""|call = ${step}
              |callInputs = ${callInputs}
              |""".stripMargin,
          minLevel = TraceLevel.VVerbose
      )
      val executableLink = getExecutableLink(step.run.name)
      // add the target step for app(let) calls
      val targetCallInput = executableLink.dxExec match {
        case _: DxWorkflow => Map.empty
        case _             => Map(Target -> (TargetParam.dxType, VString(target)))
      }
      val callInputsIR = CwlUtils.toIR(callInputs) ++ targetCallInput
      val requirementEvaluator = RequirementEvaluator(
          block.targetRequirements,
          block.targetHints,
          callInputs.map {
            case (dxName, tv) => dxName.decoded -> tv
          },
          jobMeta.workerPaths,
          step.run.inputs.map(i => i.name -> i).toMap,
          fileResolver = jobMeta.fileResolver,
          dxApi = jobMeta.dxApi
      )
      val instanceType =
        try {
          val request = requirementEvaluator.parseInstanceType
          val instanceType = jobMeta.instanceTypeDb.apply(request)
          logger.traceLimited(
              s"Precalculated instance type for ${step.run.name}: ${instanceType}"
          )
          Some(instanceType)
        } catch {
          case e: Throwable =>
            logger.traceLimited(
                s"""|Failed to precalculate the instance type for task ${step.run.name}.
                    |${e}
                    |""".stripMargin
            )
            None
        }
      val (dxExecution, execName) =
        launchJob(executableLink,
                  step.name,
                  callInputsIR,
                  nameDetail = nameDetail,
                  instanceType = instanceType.map(_.name),
                  folder = folder,
                  prefixOutputs = true)
      (dxExecution, executableLink, execName)
    }

    private def launchStepCall(stepInputs: Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))],
                               nameDetail: Option[String],
                               blockIndex: Option[Int]): (DxExecution, ExecutableLink, String) = {
      // collect all the step input values to pass to the callee
      val callInputs: Map[DxName, (CwlType, CwlValue)] = step.run.inputs.map { param =>
        val dxName = CwlDxName.fromSourceName(param.name)
        dxName -> stepInputs
          .get(dxName)
          .map(_._2)
          .orElse(Option.when(param.default.isDefined)((param.cwlType, param.default.get)))
          .orElse(Option.when(CwlOptional.isOptional(param.cwlType))((param.cwlType, NullValue)))
          .getOrElse {
            throw new Exception(
                s"missing required input ${dxName} to process ${step.run.name} at step ${step.name}"
            )
          }
      }.toMap
      launchCall(callInputs, nameDetail, folder = blockIndex.map(_.toString))
    }

    override protected def launchBlockCall(blockIndex: Int): Map[DxName, ParameterLink] = {
      assert(step.scatter.isEmpty)
      assert(step.when.isEmpty)
      val (dxExecution, executableLink, callName) =
        launchStepCall(evaluateCallStepInputs(), nameDetail = None, blockIndex = Some(blockIndex))
      jobMeta
        .createExecutionOutputLinks(dxExecution, executableLink.outputs, Some(callName))
    }

    override protected def launchConditional(): Map[DxName, ParameterLink] = {
      assert(step.scatter.isEmpty)
      assert(step.when.nonEmpty)
      val stepInputs = evaluateCallStepInputs(isConditional = true)
      val ctx = EvaluatorContext(inputs = createEvalInputs(stepInputs.values.toMap))
      val (_, cond) = eval.evaluate(step.when.get, CwlBoolean, ctx, coerce = true)
      cond match {
        case BooleanValue(true) =>
          val (dxExecution, executableLink, callName) =
            launchStepCall(stepInputs, nameDetail = None, blockIndex = Some(block.index))
          jobMeta
            .createExecutionOutputLinks(dxExecution, executableLink.outputs, Some(callName))
            .map {
              case (key, link) => key -> link.makeOptional
            }
        case _ => Map.empty
      }
    }

    private def getScatterName(items: Vector[CwlValue]): String = {
      def formatItem(item: CwlValue): String = {
        item match {
          case f: FileValue if f.location.isDefined => truncate(getFileName(f.location.get))
          case f: FileValue if f.path.isDefined     => truncate(f.path.get)
          case d: DirectoryValue if d.location.isDefined =>
            truncate(getFileName(d.location.get))
          case d: DirectoryValue if d.path.isDefined => truncate(d.path.get)
          case a: ArrayValue =>
            val itemStr =
              WorkflowExecutor.getComplexScatterName(a.items.iterator.map(i => Some(formatItem(i))))
            s"[${itemStr}]"
          case o: ObjectValue =>
            val memberStr = WorkflowExecutor.getComplexScatterName(
                o.fields.iterator.map {
                  case (k, v) => Some(s"${k}: ${formatItem(v)}")
                }
            )
            s"{${memberStr}}"
          case _ => CwlUtils.prettyFormatValue(item)
        }
      }
      items.map(formatItem).mkString(",")
    }

    override protected lazy val outputTypes: Map[DxName, Type] = block.outputs.map { param =>
      param.name -> CwlUtils.toIRType(param.cwlType)
    }.toMap

    override protected def launchScatter(): Map[DxName, ParameterLink] = {
      assert(step.scatter.nonEmpty)
      // Evaluate all the step inputs. There may be step inputs that are not passed to the callee
      // but are referred to in expressions of other step inputs' valueFrom fields. Note that the
      // step input's type here might not match the callee's input type if valueFrom is specified.
      val stepInputs = step.inputs.map(evaluateStepInput(_, cwlEnv)).toMap
      if (logger.isVerbose) {
        logger.trace(s"Workflow step ${step.name} scattering over process ${step.run.name}")
        logger.traceLimited(s"Inputs:\n  ${stepInputs.mkString("\n  ")}")
      }

      def getScatterValue(name: DxName): Option[(CwlType, Vector[CwlValue])] = {
        stepInputs.get(name) match {
          case Some((_, (t, NullValue))) if CwlOptional.isOptional(t) => None
          case Some((_, (_, ArrayValue(items)))) if items.isEmpty     => None
          case Some((_, (array: CwlArray, ArrayValue(items)))) =>
            Some(array.itemType, items)
          case Some((_, (t, ArrayValue(items))))
              if CwlOptional.unwrapOptional(t).isInstanceOf[CwlArray] =>
            Some(CwlOptional.unwrapOptional(t), items)
          case None =>
            throw new Exception(s"scatter parameter ${name} is missing from step inputs")
          case Some((_, (t, _))) =>
            throw new Exception(s"scatter parameter ${name} type ${t} is not an array")
        }
      }

      val scatterValues = scatterNames.map(getScatterValue)
      // when scatterMethod is nested_crossproduct, the output array is nested with one level of
      // depth for each scatter variable
      val outputShape = Option.when(step.scatterMethod.contains(ScatterMethod.NestedCrossproduct))(
          scatterValues.map {
            case None         => 0
            case Some((_, v)) => v.size
          }
      )
      if (scatterValues.exists(_.isEmpty)) {
        // at least one of the arrays is empty, so the scatter has no results
        return createEmptyScatterOutputs(outputShape)
      }
      val (scatterItemTypes, scatterArrays) = scatterValues.flatten.unzip
      // transform the input arrays into a single collection based on the scatter method
      val scatterCollection = step.scatterMethod match {
        case _ if scatterArrays.size == 1 => scatterArrays.transpose
        case Some(ScatterMethod.Dotproduct) if scatterArrays.map(_.size).toSet.size == 1 =>
          scatterArrays.transpose
        case Some(ScatterMethod.FlatCrossproduct) | Some(ScatterMethod.NestedCrossproduct) =>
          scatterArrays.foldLeft(Vector(Vector.empty[CwlValue])) {
            case (accu, v) => accu.flatMap(i => v.map(j => i :+ j))
          }
        case _ =>
          throw new Exception(
              "'scatterMethod' is required when there is more than one scatter variable"
          )
      }
      // select the subset of the collection for the current scatter chunk
      val (chunkCollection, next) =
        if (jobMeta.scatterStart == 0 && scatterCollection.size <= jobMeta.scatterSize) {
          (scatterCollection, None)
        } else {
          val scatterEnd = jobMeta.scatterStart + jobMeta.scatterSize
          if (scatterEnd < scatterCollection.size) {
            (scatterCollection.slice(jobMeta.scatterStart, scatterEnd), Some(scatterEnd))
          } else {
            (scatterCollection.drop(jobMeta.scatterStart), None)
          }
        }
      // the step inputs that are not listed in the 'scatter' array
      val nonScatterStepInputs = stepInputs.filterNot(i => scatterNames.contains(i._1))

      // evaluate valueFroms if necessary, get the run input values, and launch the scatter job
      def launchScatterJob(scatterValues: Vector[CwlValue]): DxExecution = {
        val allStepInputs: Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))] =
          nonScatterStepInputs ++ scatterNames
            .zip(scatterItemTypes.zip(scatterValues))
            .map {
              case (dxName, tv) => dxName -> (stepInputs(dxName)._1, tv)
            }
            .toMap
        val scatterJobInputs =
          evaluateCallStepInputs(Some(allStepInputs), isConditional = step.when.nonEmpty)
        val (dxExecution, _, _) = launchStepCall(stepInputs = scatterJobInputs,
                                                 nameDetail = Some(getScatterName(scatterValues)),
                                                 blockIndex = None)
        dxExecution
      }

      val (childJobs, skippedIndices) = if (step.when.nonEmpty) {
        // Apply `when`, which may lead to some scatter collection items being null, for which we
        // skip launching the scatter job. We keep track of which indices are null so that we can
        // reconstruct the array with the correct shape later.
        val scatterStepInputs =
          scatterNames.map(dxName => stepInputs(dxName)._1).zip(scatterItemTypes)
        val (childJobs, skippedIndices) =
          chunkCollection.foldLeft(Vector.empty[DxExecution], Vector.empty[Int]) {
            case ((execAccu, skipAccu), jobValues) =>
              val allInputs = nonScatterStepInputs.values.toMap ++ scatterStepInputs
                .zip(jobValues)
                .map {
                  case ((param, t), v) => param -> (t, v)
                }
              val ctx = EvaluatorContext(inputs = createEvalInputs(allInputs))
              val (_, cond) = eval.evaluate(step.when.get, CwlBoolean, ctx, coerce = true)
              cond match {
                case BooleanValue(true) =>
                  (execAccu :+ launchScatterJob(jobValues), skipAccu)
                case BooleanValue(false) =>
                  (execAccu, skipAccu :+ nextScatterChunkIndex)
                case other =>
                  throw new Exception(s"invalid 'when' value ${other}")
              }
          }
        (childJobs, Option.when(skippedIndices.nonEmpty)(skippedIndices))
      } else {
        (chunkCollection.map(launchScatterJob), None)
      }
      // TODO: if childJobs is empty and next is non-empty, we could just launch the next chunk
      //  from this job rather than launching a continue job
      if (childJobs.isEmpty && next.isEmpty) {
        // no jobs were launched
        createEmptyScatterOutputs(outputShape)
      } else {
        // if there are remaining chunks this calls a continue sub-job,
        // otherwise it calls a collect sub-job
        launchNextScatterChunk(childJobs, next, outputShape, skippedIndices)
      }
    }

    override protected def getScatterOutputs(
        execOutputs: Vector[Option[Map[DxName, JsValue]]],
        execName: Option[String]
    ): (Map[DxName, (Type, Value)], Option[Map[DxName, (Type, Value)]]) = {
      @tailrec
      def nestArrayTypes(itemType: Type, depth: Int): Type = {
        if (depth > 0) {
          nestArrayTypes(TArray(itemType), depth - 1)
        } else {
          itemType
        }
      }
      def nestArrayValues(array: VArray, sizes: Vector[Int]): Value = {
        if (sizes.size == 1) {
          assert(array.items.size == sizes.head)
          array
        } else {
          VArray(
              array.items
                .grouped(array.items.size / sizes.head)
                .map(group => nestArrayValues(VArray(group), sizes.tail))
                .toVector
          )
        }
      }
      val targetOutputs = step.run.outputs.map { param =>
        CwlDxName.fromSourceName(param.name) -> param.cwlType
      }.toMap
      val scatterOutputs = step.outputs.map { out =>
        val dxName = CwlDxName.fromSourceName(out.name)
        val itemType = if (step.when.nonEmpty) {
          CwlUtils.toIRType(CwlOptional.ensureOptional(targetOutputs(dxName)))
        } else {
          CwlUtils.toIRType(targetOutputs(dxName))
        }
        val flatArrayValue = createScatterOutputArray(
            execOutputs,
            CwlDxName.fromSourceName(out.name),
            itemType,
            execName
        )
        val (arrayType, arrayValue) = jobMeta.scatterOutputShape
          .map { sizes =>
            // reshape the flat array into a nested array of the specified shape
            (nestArrayTypes(itemType, sizes.size), nestArrayValues(flatArrayValue, sizes))
          }
          .getOrElse((TArray(itemType), flatArrayValue))
        dxName.pushDecodedNamespace(step.name) -> (arrayType, arrayValue)
      }.toMap
      (scatterOutputs, None)
    }
  }

  override protected def evaluateBlockInputs(
      jobInputs: Map[DxName, (Type, Value)]
  ): CwlBlockContext = {
    val block = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), jobMeta.blockPath)
    val env: Map[DxName, (CwlType, CwlValue)] = block.inputs.map { inp =>
      (inp, DxName.lookup(inp.name, jobInputs)) match {
        case (_, Some((_, (_, irValue)))) =>
          inp.name -> CwlUtils.fromIRValueWithType(irValue,
                                                   inp.cwlType,
                                                   inp.name.decoded,
                                                   isInput = true)
        case (OptionalBlockInput(name, cwlType, param), None) =>
          name -> (cwlType, NullValue)
        case _ => throw new Exception(s"missing required input ${inp.name}")
      }
    }.toMap
    CwlBlockContext(block, CwlEnv(env))
  }
}
