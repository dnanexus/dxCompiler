package dx.executor.cwl

import dx.api.DxExecution
import dx.core.Constants
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.core.ir.{Block, DxName, ExecutableLink, ParameterLink, Type, TypeSerde, Value, ValueSerde}
import dx.core.languages.Language
import dx.core.languages.cwl.{
  CwlBlock,
  CwlBlockInput,
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
  CwlOptional,
  CwlRecord,
  CwlType,
  CwlValue,
  DirectoryValue,
  Evaluator,
  EvaluatorContext,
  FileValue,
  HintUtils,
  Identifiable,
  InputParameter,
  LinkMergeMethod,
  Loadable,
  NullValue,
  ObjectValue,
  Parser,
  ParserResult,
  PickValueMethod,
  ScatterMethod,
  Workflow,
  WorkflowOutputParameter,
  WorkflowStepInput
}
import dx.executor.{JobMeta, WorkflowExecutor}
import dx.util.TraceLevel
import spray.json._

import java.net.URI

object CwlWorkflowExecutor {
  def create(jobMeta: JobMeta, separateOutputs: Boolean): CwlWorkflowExecutor = {
    // when parsing a packed workflow as a String, we need to use a baseuri -
    // it doesn't matter what it is
    val parser = Parser.create(Some(URI.create("file:/null")), hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(jobMeta.sourceCode) match {
      case Some((version, "Workflow")) if Language.parse(version) == Language.CwlV1_2 =>
        ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL Workflow document of a supported version
               |${jobMeta.sourceCode}""".stripMargin
        )
    }
    val wfName = jobMeta.getExecutableAttribute("name") match {
      case Some(JsString(name)) => name
      case _                    => throw new Exception("missing executable name")
    }
    val workflow =
      parser.parseString(jobMeta.sourceCode, defaultFrag = Some(wfName), isPacked = true) match {
        case ParserResult(wf: Workflow, _, _, _) => wf
        case other =>
          throw new Exception(s"expected CWL document to contain a Workflow, not ${other}")
      }
    CwlWorkflowExecutor(workflow, jobMeta, separateOutputs)
  }
}

case class CwlWorkflowExecutor(workflow: Workflow, jobMeta: JobMeta, separateOutputs: Boolean)
    extends WorkflowExecutor[CwlBlock](jobMeta, separateOutputs) {
  private val logger = jobMeta.logger

  override val executorName: String = "dxExecutorCwl"

  override protected lazy val typeAliases: Map[String, Type.TSchema] = {
    HintUtils.getSchemaDefs(workflow.requirements).collect {
      case (name, schema: CwlRecord) => name -> CwlUtils.toIRSchema(schema)
    }
  }

  override protected def evaluateInputs(
      jobInputs: Map[DxName, (Type, Value)]
  ): Map[DxName, (Type, Value)] = {
    // this might be the input for the entire workflow or just a subblock
    val inputs = jobMeta.blockPath match {
      case Vector() =>
        if (logger.isVerbose) {
          logger.trace(
              s"""input parameters:
                 |${workflow.inputs
                   .map { inp =>
                     s"  ${CwlUtils.prettyFormatType(inp.cwlType)} ${inp.name}"
                   }
                   .mkString("\n")}""".stripMargin
          )
        }
        jobInputs
      case path =>
        val block: CwlBlock = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), path)
        val inputTypes = block.inputs.map {
          case RequiredBlockInput(name, param) =>
            name -> CwlUtils.toIRType(param.cwlType)
          case OptionalBlockInput(name, param) =>
            val irType = CwlUtils.toIRType(param.cwlType)
            name -> Type.ensureOptional(irType)
        }.toMap
        if (logger.isVerbose) {
          logger.trace(
              s"""input parameters:
                 |${inputTypes
                   .map {
                     case (name, irType) =>
                       s"  ${TypeSerde.toString(irType)} ${name}"
                   }
                   .mkString("\n")}""".stripMargin
          )
        }
        jobInputs.collect {
          case (name, (_, v)) if inputTypes.contains(name) =>
            // coerce the input value to the target type
            val irType = inputTypes(name)
            name -> (irType, coerceTo(v, irType))
          case i => i
        }
    }
    if (logger.isVerbose) {
      logger.trace(
          s"""input values:
             |${inputs
               .map {
                 case (name, (t, v)) =>
                   s"${TypeSerde.toString(t)} ${name} = ${ValueSerde.toString(v)}}"
               }
               .mkString("\n")}""".stripMargin
      )
    }
    inputs
  }

  private case class CwlEnv[T](env: Map[DxName, T]) {
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
    logger.trace(s"Evaluating output parameters:\n  ${outputParams.mkString("\n  ")}")
    val env = CwlEnv(jobInputs)
    val irOutputs: Map[DxName, (Type, Value)] = outputParams.map {
      case (dxName, param: WorkflowOutputParameter, _) if param.sources.nonEmpty =>
        def itemsToArray(items: Vector[(Type, Value)]): (Type, Value) = {
          val (types, values) = items.unzip
          (TArray(Type.merge(types)), VArray(values))
        }
        val sourceValues = param.sources.map(src =>
          env.lookup(CwlDxName.fromDecodedName(src.frag.get)).getOrElse((TMulti.Any, VNull))
        )
        val (irType, irValue) = if (param.linkMerge.nonEmpty || param.pickValue.nonEmpty) {
          val mergedValues = if (param.linkMerge.contains(LinkMergeMethod.MergeFlattened)) {
            sourceValues.flatMap {
              case (TArray(itemType, _), VArray(items)) => items.map(i => (itemType, i))
              case (t, value)                           => Vector((t, value))
            }
          } else {
            sourceValues
          }
          if (param.pickValue.nonEmpty) {
            val nonNull = mergedValues.filterNot(_._2 == VNull)
            param.pickValue.get match {
              case PickValueMethod.FirstNonNull =>
                nonNull.headOption.getOrElse(
                    throw new Exception(
                        s"all source values are null for parameter ${param.name}"
                    )
                )
              case PickValueMethod.TheOnlyNonNull if nonNull.size == 1 => nonNull.head
              case PickValueMethod.TheOnlyNonNull =>
                throw new Exception(
                    s"there is not exactly one non-null value for parameter ${param.name}"
                )
              case PickValueMethod.AllNonNull => itemsToArray(nonNull)
            }
          } else {
            itemsToArray(mergedValues)
          }
        } else if (sourceValues.size == 1) {
          sourceValues.head
        } else {
          itemsToArray(sourceValues)
        }
        dxName -> (irType, irValue)
      case (dxName, _, cwlType) if jobInputs.contains(dxName) =>
        val irType = CwlUtils.toIRType(cwlType)
        val irValue = coerceTo(jobInputs(dxName)._2, irType)
        dxName -> (irType, irValue)
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
      CwlDxName.fromDecodedName(i.id.flatMap(_.frag).get) -> i
    }.toMap
    private lazy val eval = Evaluator.create(block.targetRequirements, block.targetHints)

    override lazy val env: Map[DxName, (Type, Value)] = CwlUtils.toIR(cwlEnv.env)

    // Each step input has zero or more sources. If there is no source, then the default value is
    // used or the value is null. Otherwise, linkMerge and pickValue are applied if specified.
    private def evaluateStepInput(
        stepInput: WorkflowStepInput,
        fullEnv: CwlEnv[(CwlType, CwlValue)]
    ): (DxName, (WorkflowStepInput, (CwlType, CwlValue))) = {
      val dxName = CwlDxName.fromSourceName(stepInput.name)
      val (cwlType, cwlValue) = if (stepInput.sources.nonEmpty) {
        def itemsToArray(items: Vector[(CwlType, CwlValue)]): (CwlType, CwlValue) = {
          val (types, values) = items.unzip
          (CwlArray(CwlType.flatten(types.distinct)), ArrayValue(values))
        }
        val sourceValues =
          stepInput.sources.map(src => fullEnv(CwlDxName.fromDecodedName(src.frag.get)))
        val (cwlType, cwlValue) =
          if (stepInput.linkMerge.nonEmpty || stepInput.pickValue.nonEmpty) {
            val mergedValues = if (stepInput.linkMerge.contains(LinkMergeMethod.MergeFlattened)) {
              sourceValues.flatMap {
                case (array: CwlArray, ArrayValue(items)) => items.map(i => (array.itemType, i))
                case (t, value)                           => Vector((t, value))
              }
            } else {
              sourceValues
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
      } else {
        (CwlAny, stepInput.default.getOrElse(NullValue))
      }
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

    private def evaluateAllStepInputs(
        extraEnv: Map[DxName, (CwlType, CwlValue)] = Map.empty
    ): Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))] = {
      val fullEnv = cwlEnv.update(extraEnv)
      // Evaluate all the step inputs. There may be step inputs that are not passed to the callee
      // but are referred to in expressions of other step inputs' valueFrom fields. Note that the
      // step input's type here might not match the callee's input type if valueFrom is specified.
      val stepInputs = step.inputs.map(evaluateStepInput(_, fullEnv)).toMap
      // evaluate valueFrom expressions - create the evaluator lazily since it might not be needed
      lazy val evalInputs: ObjectValue = createEvalInputs(stepInputs.values.toMap)
      val stepInputValues = stepInputs.map {
        case (dxName, (param, (sourceType, sourceValue))) =>
          val cwlType = runInputs.get(dxName).map(_.cwlType).getOrElse(sourceType)
          try {
            if (param.valueFrom.isDefined) {
              val ctx = EvaluatorContext(sourceValue, evalInputs)
              dxName -> (param, eval.evaluate(param.valueFrom.get, cwlType, ctx))
            } else {
              dxName -> (param, (cwlType, sourceValue.coerceTo(cwlType)))
            }
          } catch {
            case cause: Throwable =>
              throw new Exception(
                  s"""effective type ${sourceType} of input value ${sourceValue} to parameter
                     |${dxName} is not coercible to expected type ${cwlType}""".stripMargin
                    .replaceAll("\n", " "),
                  cause
              )
          }
      }
      if (logger.isVerbose) {
        logger.trace(s"Workflow step ${step.name} calling process ${step.run.name}")
        logger.traceLimited(s"Inputs:\n  ${stepInputValues.mkString("\n  ")}")
      }
      stepInputValues
    }

    private def evaluateCallInputs(
        stepInputs: Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))]
    ): Map[DxName, (CwlType, CwlValue)] = {
      // collect all the step input values to pass to the callee
      step.run.inputs.map { param =>
        val dxName = CwlDxName.fromSourceName(param.name)
        dxName -> stepInputs
          .get(dxName)
          .map(_._2)
          .orElse(Option.when(CwlOptional.isOptional(param.cwlType))((param.cwlType, NullValue)))
          .getOrElse {
            throw new Exception(
                s"missing required input ${dxName} to process ${step.run.name} at step ${step.name}"
            )
          }
      }.toMap
    }

    private def launchCall(
        callInputs: Map[DxName, (CwlType, CwlValue)],
        nameDetail: Option[String] = None,
        folder: Option[String] = None
    ): (DxExecution, ExecutableLink, String) = {
      logger.traceLimited(
          s"""|call = ${step}
              |callInputs = ${callInputs}
              |""".stripMargin,
          minLevel = TraceLevel.VVerbose
      )
      val executableLink = getExecutableLink(step.run.name)
      val targetCallInput = Map(Target -> (TargetParam.dxType, VString(block.target.name)))
      val callInputsIR = CwlUtils.toIR(callInputs) ++ targetCallInput
      val requirementEvaluator = RequirementEvaluator(
          block.targetRequirements,
          block.targetHints,
          callInputs.map {
            case (dxName, tv) => dxName.decoded -> tv
          },
          jobMeta.workerPaths,
          step.run.inputs.map(i => i.name -> i).toMap,
          dxApi = jobMeta.dxApi
      )
      val instanceType =
        try {
          val request = requirementEvaluator.parseInstanceType
          val instanceType = jobMeta.instanceTypeDb.apply(request)
          logger.traceLimited(s"Precalculated instance type for ${step.run.name}: ${instanceType}")
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
                  nameDetail,
                  instanceType.map(_.name),
                  folder = folder,
                  prefixOutputs = true)
      (dxExecution, executableLink, execName)
    }

    private def launchCall(stepInputs: Map[DxName, (WorkflowStepInput, (CwlType, CwlValue))],
                           blockIndex: Int): Map[DxName, ParameterLink] = {
      val callInputs = evaluateCallInputs(stepInputs)
      val (dxExecution, executableLink, callName) =
        launchCall(callInputs, folder = Some(blockIndex.toString))
      jobMeta.createExecutionOutputLinks(dxExecution, executableLink.outputs, Some(callName))
    }

    override protected def launchCall(blockIndex: Int): Map[DxName, ParameterLink] = {
      assert(step.scatter.isEmpty)
      assert(step.when.isEmpty)
      launchCall(evaluateAllStepInputs(), blockIndex)
    }

    override protected def launchConditional(): Map[DxName, ParameterLink] = {
      assert(step.scatter.isEmpty)
      assert(step.when.nonEmpty)
      val stepInputs = evaluateAllStepInputs()
      val ctx = EvaluatorContext(inputs = createEvalInputs(stepInputs.values.toMap))
      val (_, cond) = eval.evaluate(step.when.get, CwlBoolean, ctx)
      cond match {
        case BooleanValue(true) =>
          launchCall(stepInputs, block.index).map {
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

    override protected def launchScatter(): Map[DxName, ParameterLink] = {
      assert(step.scatter.nonEmpty)
      // construct the scatter inputs and launch one job for each, such that the inputs for each
      // job are an array of one element for each scattered variable
      val stepInputs = evaluateAllStepInputs()

      def getScatterValue(name: DxName): Option[(CwlType, Vector[CwlValue])] = {
        stepInputs.get(name) match {
          case Some((_, (t, NullValue))) if CwlOptional.isOptional(t) => None
          case Some((_, (_, ArrayValue(items)))) if items.isEmpty     => None
          case Some((_, (array: CwlArray, ArrayValue(items)))) =>
            Some(array.itemType, items)
          case None =>
            throw new Exception(s"scatter parameter ${name} is missing from env")
          case _ =>
            throw new Exception(s"scatter parameter ${name} not of type array")
        }
      }

      val scatterNames = step.scatter.map(src => CwlDxName.fromSourceName(src.name.get))
      val scatterValues = scatterNames.map(getScatterValue)
      if (scatterValues.exists(_.isEmpty)) {
        // at least one of the arrays is empty, so the scatter has no results
        return Map.empty
      }
      val (itemTypes, arrays) = scatterValues.flatten.unzip
      val scatterCollection = step.scatterMethod match {
        case ScatterMethod.Dotproduct if arrays.map(_.size).toSet.size == 1 =>
          arrays.transpose
        case ScatterMethod.FlatCrossproduct | ScatterMethod.NestedCrossproduct =>
          arrays.foldLeft(Vector(Vector.empty[CwlValue])) {
            case (accu, v) => accu.flatMap(i => v.map(j => i :+ j))
          }
      }
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

      val finalChunkCollection = if (step.when.nonEmpty) {
        val scatterInputs = scatterNames.map(dxName => stepInputs(dxName)._1).zip(itemTypes)
        val nonScatterInputs = stepInputs.keySet
          .diff(scatterNames.toSet)
          .map(dxName => stepInputs(dxName))
          .toMap
        chunkCollection.flatMap { jobValues =>
          val allInputs = nonScatterInputs ++ scatterInputs
            .zip(jobValues)
            .map {
              case ((param, t), v) => param -> (t, v)
            }
          val ctx = EvaluatorContext(inputs = createEvalInputs(allInputs))
          val (_, cond) = eval.evaluate(step.when.get, CwlBoolean, ctx)
          cond match {
            case BooleanValue(true) => Some(jobValues)
            case _                  => None
          }
        }
      } else {
        chunkCollection
      }
      val nonScatterCallInputs = stepInputs.keySet
        .diff(scatterNames.toSet)
        .map(dxName => dxName -> stepInputs(dxName)._2)
        .toMap
      val childJobs = finalChunkCollection.map { item =>
        val scatterCallInputs: Map[DxName, (CwlType, CwlValue)] = scatterNames
          .zip(itemTypes.zip(item))
          .toMap
        val callNameDetail = getScatterName(item)
        val (dxExecution, _, _) =
          launchCall(scatterCallInputs ++ nonScatterCallInputs, Some(callNameDetail))
        dxExecution
      }
      val outputShape = Option.when(step.scatterMethod.contains(ScatterMethod.NestedCrossproduct))(
          arrays.map(_.size)
      )
      next match {
        case Some(index) =>
          // there are remaining chunks - call a continue sub-job
          launchScatterContinue(childJobs, index, outputShape)
        case None =>
          // this is the last chunk - call collect sub-job to gather all the results
          launchScatterCollect(childJobs, outputShape)
      }
    }

    override protected def prepareScatterResults(
        dxSubJob: DxExecution
    ): Map[DxName, ParameterLink] = {
      // Return JBORs for all the outputs. Since the signature of the sub-job is exactly the same
      // as the parent, we can immediately exit the parent job.
      val resultTypes: Map[DxName, Type] = block.outputs.map { param =>
        param.name -> CwlUtils.toIRType(param.cwlType)
      }.toMap
      val links = jobMeta.createExecutionOutputLinks(dxSubJob, resultTypes)
      if (logger.isVerbose) {
        logger.trace(s"Scatter results for job ${dxSubJob.id}")
        logger.traceLimited(s"Types: ${resultTypes}")
        logger.traceLimited(s"Links: ${links.mkString("\n")}")
      }
      links
    }

    private def nestArrays(array: Vector[Value], sizes: Vector[Int]): Value = {
      if (sizes.size == 1) {
        assert(array.size == sizes.head)
        VArray(array)
      } else {
        VArray(
            array
              .grouped(array.size / sizes.head)
              .map(group => nestArrays(group, sizes.tail))
              .toVector
        )
      }
    }

    override protected def getScatterOutputs(
        childOutputs: Vector[Map[DxName, JsValue]],
        execName: Option[String]
    ): Map[DxName, (Type, Value)] = {
      val targetOutputs = step.run.outputs.map { param =>
        CwlDxName.fromSourceName(param.name) -> param.cwlType
      }.toMap
      val arraySizes = jobMeta.outputShape
      step.outputs.map { out =>
        val dxName = CwlDxName.fromSourceName(out.name)
        val irType = CwlUtils.toIRType(targetOutputs(dxName))
        val arrayValue =
          (createScatterOutputArray(
               childOutputs,
               CwlDxName.fromSourceName(out.name),
               irType,
               execName
           ),
           arraySizes) match {
            case (VArray(a), Some(sizes)) => nestArrays(a, sizes)
            case (a, None)                => a
            case (other, _) =>
              throw new Exception(s"invalid array value ${other}")
          }
        dxName.pushDecodedNamespace(step.name) -> (irType, arrayValue)
      }.toMap
    }
  }

  override protected def evaluateBlockInputs(
      jobInputs: Map[DxName, (Type, Value)]
  ): CwlBlockContext = {
    val block = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), jobMeta.blockPath)
    val env: Map[DxName, (CwlType, CwlValue)] = block.inputs.map {
      case inp: CwlBlockInput if jobInputs.contains(inp.name) =>
        val (_, irValue) = jobInputs(inp.name)
        inp.name -> CwlUtils.fromIRValue(irValue, inp.cwlType, inp.name.decoded, isInput = true)
      case OptionalBlockInput(name, param) =>
        name -> (param.cwlType, NullValue)
      case param =>
        throw new Exception(s"missing required input ${param.name}")
    }.toMap
    CwlBlockContext(block, CwlEnv(env))
  }
}
