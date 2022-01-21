package dx.core.languages.cwl

import dx.core.ir.{Block, BlockKind, DxName, InputKind}
import dx.cwl.{
  CommandLineTool,
  CwlAny,
  CwlArray,
  CwlMulti,
  CwlNull,
  CwlOptional,
  CwlPath,
  CwlSchema,
  CwlType,
  ExpressionTool,
  Hint,
  Identifier,
  OutputParameter,
  Parameter,
  Requirement,
  ScatterMethod,
  Workflow,
  WorkflowInputParameter,
  WorkflowStep,
  WorkflowStepInput,
  WorkflowStepOutput
}

import scala.annotation.tailrec
import scala.collection.immutable.SeqMap
import scala.util.{Success, Try}

sealed trait CwlBlockInput {
  val name: DxName
  def source: Parameter
  val kind: InputKind.InputKind

  def cwlType: CwlType = source.cwlType
}

object CwlBlockInput {
  def create(stepInput: WorkflowStepInput,
             wfId: Identifier,
             workflowSteps: Map[String, WorkflowStep],
             workflowInputs: Map[DxName, WorkflowInputParameter]): Vector[CwlBlockInput] = {
    // sources of this step input - either a workflow input (e.g. "file1") or a step output
    // (e.g. "step1/file1")
    val sources = stepInput.sources.map { src =>
      (src, CwlDxName.fromDecodedName(src.frag))
    }
    val sourceParams = sources.foldLeft(SeqMap.empty[DxName, Parameter]) {
      case (accu, (_, name)) if accu.contains(name) => accu
      case (accu, (id, name)) if id.parent.contains(wfId.frag) =>
        val wfName = name.dropNamespaces()
        accu + (wfName -> workflowInputs(wfName))
      case (accu, (id, name)) if id.parent.exists(workflowSteps.contains) =>
        val srcStep = workflowSteps(id.parent.get)
        val param = srcStep.outputs
          .collectFirst {
            case out if out.name == id.name =>
              srcStep.run.outputs
                .collectFirst {
                  case srcOut if out.name == srcOut.name => srcOut
                }
                .getOrElse(
                    throw new Exception(
                        s"""process ${srcStep.run.simpleName} does not define output parameter 
                           |${out.name}""".stripMargin.replaceAll("\n", " ")
                    )
                )
          }
          .getOrElse(
              throw new Exception(
                  s"step ${id.parent.get} does not define output parameter ${id.name}"
              )
          )
        accu + (name -> param)
      case (accu, (_, name)) if workflowInputs.contains(name) =>
        accu + (name -> workflowInputs(name))
      case (_, (_, name)) =>
        throw new Exception(s"invalid parameter source ${name}")
    }
    val hasDefault = stepInput.default.isDefined
    sourceParams.map {
      case (dxName, param) =>
        if (hasDefault || CwlOptional.isOptional(param.cwlType)) {
          OptionalBlockInput(dxName, param)
        } else {
          RequiredBlockInput(dxName, param)
        }
    }.toVector
  }
}

case class RequiredBlockInput(name: DxName, source: Parameter) extends CwlBlockInput {
  val kind: InputKind.InputKind = InputKind.Required
}

/**
  * An input that may be omitted by the caller. In that case the value will
  * be null (or None).
  */
case class OptionalBlockInput(name: DxName, source: Parameter) extends CwlBlockInput {
  val kind: InputKind.InputKind = InputKind.Optional
}

case class CwlBlockOutput(name: DxName, cwlType: CwlType, source: OutputParameter)

object CwlBlockOutput {

  /**
    * The type that results from a workflow step depends on the number of scatter sources, the
    * scatter method, the result type of the scatter target, and whether the step is conditional.
    * - dotproduct: the input arrays must be of the same length N; the output array is single
    *   dimension of length N
    * - nested_crossproduct: there is one input array for every combination of inputs; the output
    *   array is nested at depth M, where M is the number of input sources (i.e. the number of
    *   levels of scattering)
    * - flat_crossproduct: same as nested_crossproduct, but the output array is flattened to a
    *   single dimension of length N * M
    */
  def getStepOutputType(step: WorkflowStep, itemType: CwlType): CwlType = {
    val condItemType = if (step.when.isDefined) {
      CwlOptional.ensureOptional(itemType)
    } else {
      itemType
    }
    if (step.scatter.isEmpty) {
      condItemType
    } else if (step.scatter.length == 1 || step.scatterMethod.get != ScatterMethod.NestedCrossproduct) {
      CwlArray(condItemType)
    } else {
      step.scatter.indices.foldLeft(condItemType) {
        case (t, _) => CwlArray(t)
      }
    }
  }

  def create(step: WorkflowStep,
             stepOutput: WorkflowStepOutput,
             processOutput: OutputParameter): CwlBlockOutput = {
    val dxName = CwlDxName.fromSourceName(stepOutput.name, namespace = Some(step.name))
    val cwlType = getStepOutputType(step, processOutput.cwlType)
    CwlBlockOutput(dxName, cwlType, processOutput)
  }
}

/**
  * A CwlBlock is (currently) a wrapper for a single workflow step. The block inputs
  * and outputs map directly to the step inputs and outputs, which are derived from
  * the inputs and outputs of the target process, possibly augmented by the step's
  * scatter and/or conditional fields.
  * TODO: bundle ExpressionTools with the next tool or workflow step - currently
  *  they are their own step because cwltool can only execute a single workflow
  *  step in isolation.
  */
case class CwlBlock(index: Int,
                    inputs: Vector[CwlBlockInput],
                    outputs: Vector[CwlBlockOutput],
                    steps: Vector[WorkflowStep],
                    inheritedRequirements: Vector[Requirement],
                    inheritedHints: Vector[Hint])
    extends Block[CwlBlock] {
  assert(steps.nonEmpty)
  // temporary assertion
  assert(steps.size == 1)
  // the target workflow step
  val target: WorkflowStep = steps.last
  // currently empty/unused
  val prerequisites: Vector[WorkflowStep] = steps.dropRight(1)

  /**
    * Returns true if `target` can be called directly:
    * - no valueFrom
    * - at most one input source
    * - step inputs do not require down-casting to be compatible with call inputs
    */
  lazy val targetIsSimpleCall: Boolean = {
    // Whether `from` and `to` have the same native types. Assumes we've already determined that
    // `from` and `to` are compatible.
    @tailrec
    def requireDifferentNativeTypes(from: CwlType, to: CwlType): Boolean = {
      (from, to) match {
        case (CwlNull, CwlNull)               => false
        case (CwlNull, CwlOptional(_))        => false
        case (CwlOptional(f), CwlOptional(t)) => requireDifferentNativeTypes(f, t)
        case (f, t: CwlArray)                 =>
          // it's allowed to pass P to array:P
          requireDifferentNativeTypes(f, t.itemType)
        case (CwlAny | _: CwlMulti | _: CwlPath | _: CwlSchema,
              CwlAny | _: CwlMulti | _: CwlPath | _: CwlSchema) =>
          // hash to hash
          false
        case (_, _: CwlMulti | CwlAny) => true
        case _                         => false
      }
    }

    target.run match {
      case _: ExpressionTool | _: CommandLineTool if CwlUtils.isSimpleCall(target) =>
        val blockInputs = inputs.map(i => i.name -> i).toMap
        val stepInputs = target.inputs.map(i => i.name -> i).toMap
        !target.run.inputs.exists { targetInput =>
          // Check that the target input type is compatible with the block input type.
          // The call is not simple if a downcast is required from the block input to
          // the target input (e.g. `Any` -> `File`), or if the target input must be given
          // as a hash and the block input is not (e.g. `String` -> `Any`).
          stepInputs.get(targetInput.name).exists { stepInput =>
            stepInput.sources.headOption.exists { src =>
              DxName.lookup(CwlDxName.fromDecodedName(src.frag), blockInputs).exists {
                case (_, blockInput) =>
                  val blockInputType =
                    if (stepInput.linkMerge.isEmpty || CwlUtils.isArray(blockInput.cwlType)) {
                      blockInput.cwlType
                    } else {
                      CwlArray(blockInput.cwlType)
                    }
                  val targetInputType =
                    if (stepInput.default.isDefined || targetInput.default.isDefined) {
                      CwlOptional.ensureOptional(targetInput.cwlType)
                    } else {
                      targetInput.cwlType
                    }
                  Try(CwlUtils.requiresDowncast(blockInputType, targetInputType)) match {
                    case Success(false) =>
                      requireDifferentNativeTypes(blockInputType, targetInputType)
                    case _ =>
                      // this may fail due to the need to cast an optional step input to a
                      // non-optional process input, in which case we need a fragment to perform
                      // the type casting
                      true
                  }
              }
            }
          }
        }
      case _ => false
    }
  }

  /**
    * The kind of block this is.
    */
  override lazy val kind: BlockKind.BlockKind = {
    if (target.scatter.nonEmpty) {
      BlockKind.ScatterOneCall
    } else if (target.when.nonEmpty) {
      BlockKind.ConditionalOneCall
    } else if (targetIsSimpleCall) {
      BlockKind.CallDirect
    } else {
      BlockKind.CallFragment
    }
  }

  override lazy val getName: Option[String] = {
    kind match {
      case BlockKind.CallDirect | BlockKind.CallFragment =>
        Some(s"frag ${steps.last.name}")
      case BlockKind.ScatterOneCall =>
        Some(s"scatter (${steps.last.scatter.mkString(",")})")
      case BlockKind.ConditionalOneCall =>
        Some(s"if (${CwlUtils.prettyFormatValue(steps.last.when.get)})")
      case _ => None
    }
  }

  def targetRequirements: Vector[Requirement] = {
    inheritedRequirements ++ target.requirements ++ target.run.requirements
  }

  def targetHints: Vector[Hint] = {
    inheritedHints ++ target.hints ++ target.run.hints
  }

  override def getSubBlock(index: Int): CwlBlock = {
    kind match {
      case BlockKind.ConditionalOneCall | BlockKind.ScatterOneCall =>
        target.run match {
          case _: CommandLineTool if index == 0 =>
            CwlBlock(0, inputs, outputs, Vector(target), inheritedRequirements, inheritedHints)
          case wf: Workflow =>
            val innerBlocks =
              CwlBlock.createBlocks(wf,
                                    inheritedRequirements ++ target.requirements,
                                    inheritedHints ++ target.hints)
            innerBlocks(index)
          case other =>
            throw new Exception(s"unexpected process ${other}")
        }
      case _ =>
        throw new Exception(s"block ${this} does not have inner elements")
    }
  }

  override lazy val inputNames: Set[DxName] = inputs.map(_.name).toSet

  override lazy val outputNames: Set[DxName] = outputs.map(_.name).toSet

  private def prettyFormatStep(step: WorkflowStep): String = {
    val parts = Vector(
        if (step.scatter.nonEmpty) {
          Some(s"scatter (${step.scatter.mkString(",")}) {", "}")
        } else None,
        if (step.when.nonEmpty) {
          Some(s"if (${CwlUtils.prettyFormatValue(step.when.get)}) {", "}")
        } else None,
        Some(
            step.run match {
              case tool: CommandLineTool => s"tool ${tool.simpleName}"
              case expr: ExpressionTool  => s"expr ${expr.simpleName}"
              case wf: Workflow          => s"workflow ${wf.simpleName}"
              case other =>
                throw new Exception(s"unexpected process ${other}")
            },
            ""
        )
    ).flatten
    val ((prefix, suffix), _) = parts.zipWithIndex.reduce { (x, y) =>
      val ((p1, s1), _) = x
      val ((p2, s2), i) = y
      val indent = " " * 4 * i
      val prefix = s"${p1}\n${indent}${p2}"
      val suffix = if (s2.nonEmpty) {
        s"${indent}${s2}\n${s1}"
      } else {
        s1
      }
      ((prefix, suffix), 0)
    }
    if (suffix.nonEmpty) {
      s"${prefix}\n${suffix}"
    } else {
      prefix
    }
  }

  override lazy val prettyFormat: String = {
    val inputStr = inputs.map { param =>
      s"${CwlUtils.prettyFormatType(param.cwlType)} ${param.name}"
    }
    val outputStr = outputs.map { param =>
      s"${CwlUtils.prettyFormatType(param.cwlType)} ${param.name}"
    }
    val bodyStr = steps.map(prettyFormatStep).mkString("\n")
    s"""Block(${index}, ${kind})
       |  Inputs: ${inputStr}
       |  Outputs: ${outputStr}
       |  Body: ${bodyStr}""".stripMargin
  }
}

object CwlBlock {
  def createBlocks(wf: Workflow,
                   inheritedRequirements: Vector[Requirement] = Vector.empty,
                   inheritedHints: Vector[Hint] = Vector.empty): Vector[CwlBlock] = {
    val wfInputs: Map[DxName, WorkflowInputParameter] = wf.inputs.collect {
      case i if i.hasName => CwlDxName.fromSourceName(i.name) -> i
    }.toMap

    // sort steps in dependency order
    def orderSteps(
        wfId: Identifier,
        steps: Vector[WorkflowStep],
        deps: SeqMap[String, WorkflowStep] = SeqMap.empty
    ): SeqMap[String, WorkflowStep] = {
      val (satisfied, unsatisfied) =
        steps.partition { step =>
          step.inputs.forall { inp =>
            inp.sources.forall {
              case id if id.parent.contains(wfId.frag)   => true
              case id if id.parent.exists(deps.contains) => true
              case id =>
                wfInputs.contains(CwlDxName.fromSourceName(id.name))
            }
          }
        }
      if (satisfied.isEmpty) {
        throw new Exception("could not satisfy any dependencies")
      }
      val newDeps = deps ++ satisfied.map(step => step.frag -> step).to(SeqMap)
      if (unsatisfied.isEmpty) {
        newDeps
      } else {
        newDeps ++ orderSteps(wfId, unsatisfied, newDeps)
      }
    }

    if (wf.steps.isEmpty) {
      Vector.empty[CwlBlock]
    } else {
      val orderedSteps = orderSteps(wf.id.get, wf.steps)
      orderedSteps.values.zipWithIndex.map {
        case (step, index) =>
          val blockInputs: Vector[CwlBlockInput] = {
            step.inputs
              .flatMap { inp =>
                CwlBlockInput.create(inp, wf.id.get, orderedSteps, wfInputs)
              }
              .distinct
              .foldLeft(
                  SeqMap.empty[DxName, (Option[RequiredBlockInput], Option[OptionalBlockInput])]
              )({
                case (accu, req: RequiredBlockInput) if accu.contains(req.name) =>
                  val (existingReq, existingOpt) = accu(req.name)
                  if (existingReq.nonEmpty) {
                    throw new Exception(s"multiple different required inputs for name ${req.name}")
                  }
                  accu + (req.name -> (Some(req), existingOpt))
                case (accu, opt: OptionalBlockInput) if accu.contains(opt.name) =>
                  val (existingReq, existingOpt) = accu(opt.name)
                  if (existingOpt.nonEmpty) {
                    throw new Exception(s"multiple different optional inputs for name ${opt.name}")
                  }
                  accu + (opt.name -> (existingReq, Some(opt)))
                case (accu, req: RequiredBlockInput) =>
                  accu + (req.name -> (Some(req), Option.empty[OptionalBlockInput]))
                case (accu, opt: OptionalBlockInput) =>
                  accu + (opt.name -> (Option.empty[RequiredBlockInput], Some(opt)))
                case (_, other) =>
                  throw new Exception(s"unexpected block input ${other}")
              })
              .values
              .map {
                case (req, opt) => req.orElse(opt).get
              }
              .toVector
          }
          val procOutputs = step.run.outputs.collect {
            case o if o.hasName => o.name -> o
          }.toMap
          val blockOutputs = step.outputs
            .foldLeft(SeqMap.empty[DxName, CwlBlockOutput]) {
              case (_, out) if !procOutputs.contains(out.name) =>
                throw new Exception(
                    s"step output ${out} does not match any output of callable ${step.run.id}"
                )
              case (accu, out) =>
                val blockOutput = CwlBlockOutput.create(step, out, procOutputs(out.name))
                if (accu.contains(blockOutput.name)) {
                  throw new Exception(s"duplicate output parameter name ${blockOutput.name}")
                }
                accu + (blockOutput.name -> blockOutput)
            }
            .values
            .toVector
          CwlBlock(index,
                   blockInputs,
                   blockOutputs,
                   Vector(step),
                   inheritedRequirements ++ wf.requirements,
                   inheritedHints ++ wf.hints)
      }.toVector
    }
  }
}
