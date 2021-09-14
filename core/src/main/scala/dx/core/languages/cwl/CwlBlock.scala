package dx.core.languages.cwl

import dx.core.ir.{Block, BlockKind, DxName, InputKind}
import dx.cwl.{
  CommandLineTool,
  CwlArray,
  CwlOptional,
  CwlType,
  ExpressionTool,
  Hint,
  OutputParameter,
  Parameter,
  Requirement,
  Workflow,
  WorkflowInputParameter,
  WorkflowStep,
  WorkflowStepOutput
}

import scala.collection.immutable.SeqMap

sealed trait CwlBlockInput {
  val name: DxName
  def source: Parameter
  val kind: InputKind.InputKind

  def cwlType: CwlType = source.cwlType
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
  def create(step: WorkflowStep,
             stepOutput: WorkflowStepOutput,
             targetParam: OutputParameter): CwlBlockOutput = {
    val dxName = CwlDxName.fromSourceName(stepOutput.name, namespace = Some(step.name))
    val cwlType = (step.when.isDefined, step.scatter.nonEmpty) match {
      case (true, true)   => CwlArray(CwlOptional.ensureOptional(targetParam.cwlType))
      case (true, false)  => CwlOptional.ensureOptional(targetParam.cwlType)
      case (false, true)  => CwlArray(targetParam.cwlType)
      case (false, false) => targetParam.cwlType
    }
    CwlBlockOutput(dxName, cwlType, targetParam)
  }
}

case class CwlBlock(index: Int,
                    inputs: Vector[CwlBlockInput],
                    outputs: Vector[CwlBlockOutput],
                    steps: Vector[WorkflowStep],
                    inheritedRequirements: Vector[Requirement],
                    inheritedHints: Vector[Hint])
    extends Block[CwlBlock] {
  assert(steps.nonEmpty)

  val prerequisites: Vector[WorkflowStep] = steps.dropRight(1)
  val target: WorkflowStep = steps.last

  /**
    * Returns true if `target` can be called directly:
    * - no valueFrom
    * - at most one input source
    * - step inputs do not require down-casting to be compatible with call inputs
    */
  lazy val targetIsSimpleCall: Boolean = {
    target.run match {
      case _: ExpressionTool | _: CommandLineTool if CwlUtils.isSimpleCall(target) =>
        val blockInputs = inputs.map(i => i.name -> i).toMap
        val stepInputs = target.inputs.map(i => i.name -> i).toMap
        !target.run.inputs.exists { callInput =>
          // check that the call input type is compatible with the block input type
          // the call is not simple if a downcast is required (e.g. Any -> File)
          stepInputs.get(callInput.name).exists { stepInput =>
            stepInput.sources.headOption.exists { src =>
              blockInputs.get(CwlDxName.fromDecodedName(src.frag.get)).exists { blockInput =>
                try {
                  CwlUtils.requiresDowncast(blockInput.cwlType, callInput.cwlType)
                } catch {
                  case cause: Throwable =>
                    throw new Exception(
                        s"block input ${blockInput} is not compatible with call input ${callInput}",
                        cause
                    )
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
              case tool: CommandLineTool => s"tool ${tool.name}"
              case expr: ExpressionTool  => s"expr ${expr.name}"
              case wf: Workflow          => s"workflow ${wf.name}"
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
    // TODO: bundle ExpressionTools with the next tool or workflow step - currently
    //  they are their own step because cwltool can only execute a single workflow
    //  step in isolation

    val wfInputs: Map[DxName, WorkflowInputParameter] = wf.inputs.collect {
      case i if i.hasName =>
        CwlDxName.fromSourceName(i.name) -> i
    }.toMap

    // sort steps in dependency order
    def orderSteps(
        steps: Vector[WorkflowStep],
        deps: SeqMap[String, WorkflowStep] = SeqMap.empty
    ): SeqMap[String, WorkflowStep] = {
      val (satisfied, unsatisfied) =
        steps.partition { step =>
          step.inputs.forall { inp =>
            inp.sources.forall {
              case id if id.parent.isDefined => deps.contains(id.parent.get)
              case id =>
                wfInputs.contains(CwlDxName.fromSourceName(id.name.get))
            }
          }
        }
      if (satisfied.isEmpty) {
        throw new Exception("could not satisfy any dependencies")
      }
      val newDeps = deps ++ satisfied.map(step => step.name -> step).to(SeqMap)
      if (unsatisfied.isEmpty) {
        newDeps
      } else {
        newDeps ++ orderSteps(unsatisfied, newDeps)
      }
    }

    if (wf.steps.isEmpty) {
      Vector.empty[CwlBlock]
    } else {
      val orderedSteps = orderSteps(wf.steps)
      orderedSteps.values.zipWithIndex.map {
        case (step, index) =>
          val blockInputs: Vector[CwlBlockInput] = {
            val inputSources: Vector[CwlBlockInput] = step.inputs.flatMap { inp =>
              // sources of this step input - either a workflow input
              // like "file1" or a step output like "step1/file1"
              val sources = inp.sources.map { src =>
                (src, CwlDxName.fromDecodedName(src.frag.get))
              }
              val sourceParams = sources.foldLeft(SeqMap.empty[DxName, Parameter]) {
                case (accu, (_, name)) if accu.contains(name) => accu
                case (accu, (id, name)) if id.parent.isDefined =>
                  val srcStep = orderedSteps(id.parent.get)
                  val param = srcStep.outputs
                    .collectFirst {
                      case out if out.name == id.name.get =>
                        srcStep.run.outputs
                          .collectFirst {
                            case srcOut if out.name == srcOut.name => srcOut
                          }
                          .getOrElse(
                              throw new Exception(
                                  s"""process ${srcStep.run.name} does not define output parameter 
                                     |${out.name}""".stripMargin.replaceAll("\n", " ")
                              )
                          )
                    }
                    .getOrElse(
                        throw new Exception(
                            s"step ${id.parent.get} does not define output parameter ${id.name.get}"
                        )
                    )
                  accu + (name -> param)
                case (accu, (_, name)) if wfInputs.contains(name) =>
                  accu + (name -> wfInputs(name))
                case (_, (_, name)) =>
                  throw new Exception(s"invalid parameter source ${name}")
              }
              val hasDefault = inp.default.isDefined
              sourceParams.map {
                case (dxName, param) =>
                  if (hasDefault || CwlOptional.isOptional(param.cwlType)) {
                    OptionalBlockInput(dxName, param)
                  } else {
                    RequiredBlockInput(dxName, param)
                  }
              }.toVector
            }
            inputSources.distinct
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
