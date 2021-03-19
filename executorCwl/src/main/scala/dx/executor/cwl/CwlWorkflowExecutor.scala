package dx.executor.cwl

import dx.api.DxExecution
import dx.core.Constants
import dx.core.ir.Type.TString
import dx.core.ir.Value.{VNull, VString}
import dx.core.ir.{
  Block,
  BlockKind,
  ExecutableLink,
  ParameterLink,
  Type,
  TypeSerde,
  Value,
  ValueSerde
}
import dx.core.languages.Language
import dx.core.languages.cwl.{
  CwlBlock,
  CwlBlockInput,
  CwlUtils,
  DxHintSchema,
  OptionalBlockInput,
  RequiredBlockInput,
  RequirementEvaluator
}
import dx.cwl.{
  ArrayValue,
  CwlArray,
  CwlOptional,
  CwlRecord,
  CwlType,
  CwlValue,
  DirectoryValue,
  FileValue,
  HintUtils,
  NullValue,
  ObjectValue,
  Parser,
  ScatterMethod,
  Workflow
}
import dx.executor.{JobMeta, WorkflowExecutor}
import dx.util.TraceLevel
import spray.json._

object CwlWorkflowExecutor {
  def create(jobMeta: JobMeta, separateOutputs: Boolean): CwlWorkflowExecutor = {
    val parser = Parser.create(hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(jobMeta.sourceCode) match {
      case Some((version, "CommandLineTool")) if Language.parse(version) == Language.CwlV1_2 => ()
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
    val workflow = parser.parseString(jobMeta.sourceCode, name = Some(wfName)) match {
      case tool: Workflow => tool
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
      jobInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    // This might be the input for the entire workflow or just a subblock.
    // If it is for a sublock, it may be for the body of a conditional or
    // scatter, in which case we only need the inputs of the body statements.
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
          case RequiredBlockInput(name, cwlType) =>
            name -> CwlUtils.toIRType(cwlType)
          case OptionalBlockInput(name, cwlType) =>
            val irType = CwlUtils.toIRType(cwlType)
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
            name -> (irType, Value.coerceTo(v, irType))
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

  override protected def evaluateOutputs(jobInputs: Map[String, (Type, Value)],
                                         addReorgStatus: Boolean): Map[String, (Type, Value)] = {
    // This might be the output for the entire workflow or just a subblock.
    // If it is for a sublock, it may be for the body of a conditional or
    // scatter, in which case we only need the outputs of the body statements.
    val (outputParams, optional, array) = jobMeta.blockPath match {
      case Vector() => (workflow.outputs, false, false)
      case path =>
        val block = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), path)
        block.kind match {
          case BlockKind.CallDirect | BlockKind.CallFragment =>
            (block.outputs, false, false)
          case _ =>
            val step = block.target.get
            (block.outputs, step.when.isDefined, step.scatter.nonEmpty)
        }
    }
    val irOutputs = outputParams.map {
      case param if jobInputs.contains(param.name) =>
        val paramIrType = CwlUtils.toIRType(param.cwlType)
        val irType = (optional, array) match {
          case (true, true)   => Type.TArray(Type.ensureOptional(paramIrType))
          case (true, false)  => Type.ensureOptional(paramIrType)
          case (false, true)  => Type.TArray(paramIrType)
          case (false, false) => paramIrType
        }
        param.name -> (irType, Value.coerceTo(jobInputs(param.name)._2, irType))
      case param if CwlOptional.isOptional(param.cwlType) =>
        param.name -> (CwlUtils.toIRType(param.cwlType), VNull)
      case param =>
        throw new Exception(s"missing required output ${param.name}")
    }.toMap
    if (addReorgStatus) {
      irOutputs + (Constants.ReorgStatus -> (TString, VString(Constants.ReorgStatusCompleted)))
    } else {
      irOutputs
    }
  }

  case class CwlBlockContext(block: CwlBlock, cwlEnv: Map[String, (CwlType, CwlValue)])
      extends BlockContext {
    private val step = block.target.get

    override lazy val env: Map[String, (Type, Value)] = CwlUtils.toIR(cwlEnv)

    private def evaluateCallInputs(
        extraEnv: Map[String, (CwlType, CwlValue)] = Map.empty
    ): Map[String, (CwlType, CwlValue)] = {
      val fullEnv = cwlEnv ++ extraEnv
      val stepInputs = step.inputs.collect(inp => inp.name -> inp).toMap
      step.run.inputs.flatMap {
        case param if stepInputs.contains(param.name) =>
          val step = stepInputs(param.name)
          step.source.map { src =>
            val name = src.name.get
            name -> fullEnv(name)
          }
        case param if CwlOptional.isOptional(param.cwlType) =>
          logger.trace(s"no input for optional input ${param.name} to step ${step.name}")
          Vector.empty
        case param =>
          throw new Exception(s"missing non-optional input ${param.name} to step ${step.name}")
      }.toMap
    }

    private def launchCall(
        callInputs: Map[String, (CwlType, CwlValue)],
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
      val callInputsIR = CwlUtils.toIR(callInputs)
      val requirementEvaluator = RequirementEvaluator(
          block.targetRequirements,
          block.targetHints,
          callInputs,
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

    override protected def launchCall(blockIndex: Int): Map[String, ParameterLink] = {
      val callInputs = evaluateCallInputs()
      val (dxExecution, executableLink, callName) =
        launchCall(callInputs, folder = Some(blockIndex.toString))
      jobMeta.createExecutionOutputLinks(dxExecution, executableLink.outputs, Some(callName))
    }

    override protected def launchConditional(): Map[String, ParameterLink] = {
      launchCall(block.index).map {
        case (key, link) => key -> link.makeOptional
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

    private def launchScatterCallJobs(scatterParams: Vector[String],
                                      itemTypes: Vector[CwlType],
                                      collection: Vector[Vector[CwlValue]]): Vector[DxExecution] = {
      val otherInputs = cwlEnv.view.filterKeys(!scatterParams.contains(_))
      collection.map { item =>
        val callInputs = scatterParams.zip(itemTypes.zip(item)).toMap ++ otherInputs
        val callNameDetail = getScatterName(item)
        val (dxExecution, _, _) = launchCall(callInputs, Some(callNameDetail))
        dxExecution
      }
    }

    override protected def prepareScatterResults(
        dxSubJob: DxExecution
    ): Map[String, ParameterLink] = {
      val resultTypes = block.outputs.map { param =>
        param.name -> CwlUtils.toIRType(param.cwlType)
      }.toMap
      // Return JBORs for all the outputs. Since the signature of the sub-job
      // is exactly the same as the parent, we can immediately exit the parent job.
      val links = jobMeta.createExecutionOutputLinks(dxSubJob, resultTypes)
      if (logger.isVerbose) {
        val linkStr = links.mkString("\n")
        logger.traceLimited(s"resultTypes=${resultTypes}")
        logger.traceLimited(s"promises=${linkStr}")
      }
      links
    }

    private def getScatterValues: Vector[Option[(CwlType, Vector[CwlValue])]] = {
      step.scatter.map(name =>
        cwlEnv.get(name) match {
          case Some((t, NullValue)) if CwlOptional.isOptional(t) => None
          case Some((_, ArrayValue(items))) if items.isEmpty     => None
          case Some((array: CwlArray, ArrayValue(items))) =>
            Some(array.itemType, items)
          case None =>
            throw new Exception(s"scatter parameter ${name} is missing from env")
          case _ =>
            throw new Exception(s"scatter parameter ${name} not of type array")
        }
      )
    }

    override protected def launchScatter(): Map[String, ParameterLink] = {
      // construct the scatter inputs and launch one job for each, such
      // that the inputs for each job are an array of one element for
      // each scattered variable
      val scatterValues = getScatterValues
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
      val childJobs = launchScatterCallJobs(step.scatter, itemTypes, chunkCollection)
      next match {
        case Some(index) =>
          // there are remaining chunks - call a continue sub-job
          launchScatterContinue(childJobs, index)
        case None =>
          // this is the last chunk - call collect sub-job to gather all the results
          launchScatterCollect(childJobs)
      }
    }

    private def nestArrays(array: Vector[Value], sizes: Vector[Int]): Value = {
      if (sizes.size == 1) {
        assert(array.size == sizes.head)
        Value.VArray(array)
      } else {
        Value.VArray(
            array
              .grouped(array.size / sizes.head)
              .map(group => nestArrays(group, sizes.tail))
              .toVector
        )
      }
    }

    override protected def getScatterOutputs(
        childOutputs: Vector[Map[String, JsValue]],
        execName: Option[String]
    ): Map[String, (Type, Value)] = {
      val targetOutputs = step.run.outputs.map(param => param.name -> param.cwlType).toMap
      val arraySizes = if (step.scatterMethod.get == ScatterMethod.NestedCrossproduct) {
        val (_, arrays) = getScatterValues.flatten.unzip
        Some(arrays.map(_.size))
      } else {
        None
      }
      step.outputs.map { out =>
        val irType = CwlUtils.toIRType(targetOutputs(out.name))
        val arrayValue =
          (createScatterOutputArray(childOutputs, out.name, irType, execName), arraySizes) match {
            case (Value.VArray(a), Some(sizes)) => nestArrays(a, sizes)
            case (a, None)                      => a
            case (other, _) =>
              throw new Exception(s"invalid array value ${other}")
          }
        s"${step.name}/${out.name}" -> (irType, arrayValue)
      }.toMap
    }
  }

  override protected def evaluateBlockInputs(
      jobInputs: Map[String, (Type, Value)]
  ): CwlBlockContext = {
    val block = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), jobMeta.blockPath)
    val env = block.inputs.map {
      case param: CwlBlockInput if jobInputs.contains(param.name) =>
        val (_, irValue) = jobInputs(param.name)
        param.name -> CwlUtils.fromIRValue(irValue, param.cwlType, param.name, isInput = true)
      case OptionalBlockInput(name, cwlType) =>
        name -> (cwlType, NullValue)
      case param =>
        throw new Exception(s"missing required input ${param.name}")
    }.toMap
    CwlBlockContext(block, env)
  }
}
