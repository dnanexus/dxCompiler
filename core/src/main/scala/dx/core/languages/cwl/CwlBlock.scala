package dx.core.languages.cwl

import dx.core.ir.{Block, BlockKind, DxName, InputKind}
import dx.cwl.{
  CommandLineTool,
  CwlOptional,
  CwlType,
  ExpressionTool,
  Hint,
  OutputParameter,
  Parameter,
  Requirement,
  Workflow,
  WorkflowInputParameter,
  WorkflowStep
}

import scala.collection.immutable.TreeSeqMap

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

case class CwlBlock(index: Int,
                    inputs: Map[DxName, CwlBlockInput],
                    outputs: Map[DxName, OutputParameter],
                    steps: Vector[WorkflowStep],
                    inheritedRequirements: Vector[Requirement],
                    inheritedHints: Vector[Hint])
    extends Block[CwlBlock] {
  assert(steps.nonEmpty)

  /**
    * The kind of block this is.
    */
  override lazy val kind: BlockKind.BlockKind = {
    val call = steps.last
    if (call.scatter.nonEmpty) {
      BlockKind.ScatterOneCall
    } else if (call.when.nonEmpty) {
      BlockKind.ConditionalOneCall
    } else {
      call.run match {
        case _: ExpressionTool | _: CommandLineTool if CwlUtils.isSimpleCall(call) =>
          BlockKind.CallDirect
        case _ => BlockKind.CallFragment
      }
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

  lazy val (prerequisites, target): (Vector[WorkflowStep], Option[WorkflowStep]) = {
    (steps.dropRight(1), Some(steps.last))
  }

  def targetRequirements: Vector[Requirement] = {
    target
      .map(step => inheritedRequirements ++ step.requirements ++ step.run.requirements)
      .getOrElse(Vector.empty)
  }

  def targetHints: Vector[Hint] = {
    target
      .map(step => inheritedHints ++ step.hints ++ step.run.hints)
      .getOrElse(Vector.empty)
  }

  override def getSubBlock(index: Int): CwlBlock = {
    (kind, target) match {
      case (BlockKind.ConditionalOneCall | BlockKind.ScatterOneCall, Some(step)) =>
        step.run match {
          case _: CommandLineTool if index == 0 =>
            CwlBlock(0, inputs, outputs, Vector(step), inheritedRequirements, inheritedHints)
          case wf: Workflow =>
            val innerBlocks =
              CwlBlock.createBlocks(wf,
                                    inheritedRequirements ++ step.requirements,
                                    inheritedHints ++ step.hints)
            innerBlocks(index)
          case other =>
            throw new Exception(s"unexpected process ${other}")
        }
      case _ =>
        throw new Exception(s"block ${this} does not have inner elements")
    }
  }

  override lazy val inputNames: Set[DxName] = inputs.keySet

  override lazy val outputNames: Set[DxName] = outputs.keySet

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
    val inputStr = inputs.values.map { param =>
      s"${CwlUtils.prettyFormatType(param.cwlType)} ${param.name}"
    }
    val outputStr = outputs.values.map { param =>
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
        deps: TreeSeqMap[String, WorkflowStep] = TreeSeqMap.empty
    ): TreeSeqMap[String, WorkflowStep] = {
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
      val newDeps = deps ++ satisfied.map(step => step.name -> step).to(TreeSeqMap)
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
          val blockInputs: Map[DxName, CwlBlockInput] = {
            val inputSources: Vector[Map[DxName, CwlBlockInput]] = step.inputs.map { inp =>
              // sources of this step input - either a workflow input
              // like "file1" or a step output like "step1/file1"
              val sources = inp.sources.map { src =>
                (src, CwlDxName.fromDecodedName(src.frag.get))
              }
              val sourceParams = sources.foldLeft(Map.empty[DxName, Parameter]) {
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
                                  s"process ${srcStep.run.name} does not define output parameter ${out.name}"
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
                    dxName -> OptionalBlockInput(dxName, param)
                  } else {
                    dxName -> RequiredBlockInput(dxName, param)
                  }
              }
            }

            inputSources.flatMap(_.toVector).groupBy(_._1).map {
              case (name, inputs) =>
                val param = inputs.map(_._2).distinct match {
                  case Vector(i)                                              => i
                  case Vector(req: RequiredBlockInput, _: OptionalBlockInput) => req
                  case Vector(_: OptionalBlockInput, req: RequiredBlockInput) => req
                  case _ =>
                    throw new Exception(
                        s"multiple incompatible inputs ${inputs} map to name ${name}"
                    )
                }
                name -> param
            }
          }
          val procOutputs = step.run.outputs.collect {
            case o if o.hasName => o.name -> o
          }.toMap
          val blockOutputs = step.outputs
            .map {
              case out if procOutputs.contains(out.name) =>
                val dxName = CwlDxName.fromSourceName(out.name, namespace = Some(step.name))
                dxName -> out
              case out =>
                throw new Exception(s"invalid output parameter name ${out.name}")
            }
            .foldLeft(Map.empty[DxName, OutputParameter]) {
              case (accu, (dxName, out)) if !accu.contains(dxName) =>
                accu + (dxName -> procOutputs(out.name))
              case (_, (_, out)) =>
                throw new Exception(s"duplicate output parameter name ${out.name}")
            }
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
