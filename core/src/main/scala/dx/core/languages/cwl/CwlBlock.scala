package dx.core.languages.cwl

import dx.core.ir.{Block, BlockKind, InputKind, Parameter => IRParameter}
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
  WorkflowStep
}

import scala.collection.immutable.TreeSeqMap

sealed trait CwlBlockInput {
  val name: String
  val cwlType: CwlType
  val kind: InputKind.InputKind
}

case class RequiredBlockInput(name: String, cwlType: CwlType) extends CwlBlockInput {
  val kind: InputKind.InputKind = InputKind.Required
}

/**
  * An input that may be omitted by the caller. In that case the value will
  * be null (or None).
  */
case class OptionalBlockInput(name: String, cwlType: CwlType) extends CwlBlockInput {
  val kind: InputKind.InputKind = InputKind.Optional
}

case class CwlBlock(index: Int,
                    inputs: Map[String, CwlBlockInput],
                    outputs: Map[String, OutputParameter],
                    steps: Vector[WorkflowStep],
                    inheritedRequirements: Vector[Requirement],
                    inheritedHints: Vector[Hint])
    extends Block[CwlBlock] {
  assert(steps.nonEmpty)

  /**
    * The kind of block this is.
    */
  override lazy val kind: BlockKind.BlockKind = {
    (steps.last, steps.last.run) match {
      case (_, _: ExpressionTool) =>
        BlockKind.ExpressionsOnly
      case (step, _: CommandLineTool) if step.scatter.isEmpty && step.when.isEmpty =>
        BlockKind.CallDirect
      case (step, _) if step.scatter.nonEmpty =>
        BlockKind.ScatterOneCall
      case (step, _) if step.when.nonEmpty =>
        BlockKind.ConditionalOneCall
      case _ =>
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

  lazy val (prerequisites, target): (Vector[WorkflowStep], Option[WorkflowStep]) = {
    kind match {
      case BlockKind.ExpressionsOnly => (steps, None)
      case _                         => (steps.dropRight(1), Some(steps.last))
    }
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

  override lazy val inputNames: Set[String] = inputs.keySet

  override lazy val outputNames: Set[String] = outputs.values.map(_.name).toSet

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
    val inputStr = inputs.map {
      case (name, param) =>
        s"${CwlUtils.prettyFormatType(param.cwlType)} ${name}"
    }
    val outputStr = outputs.map {
      case (name, param) =>
        s"${CwlUtils.prettyFormatType(param.cwlType)} ${name}"
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

    val wfInputs = wf.inputs.collect {
      case i if i.hasName => i.name -> i
    }.toMap

    // sort steps in dependency order
    def orderSteps(
        steps: Vector[WorkflowStep],
        deps: TreeSeqMap[String, WorkflowStep] = TreeSeqMap.empty
    ): TreeSeqMap[String, WorkflowStep] = {
      val (satisfied, unsatisfied) =
        steps.partition { step =>
          step.inputs.forall { inp =>
            inp.source.forall {
              case id if id.parent.isDefined => deps.contains(id.parent.get)
              case id                        => wfInputs.contains(id.name.get)
            }
          }
        }
      val newDeps = deps ++ satisfied.map(step => step.name -> step).to(TreeSeqMap)
      newDeps ++ orderSteps(unsatisfied, newDeps)
    }

    val orderedSteps = orderSteps(wf.steps)

    orderedSteps.values.zipWithIndex.map {
      case (step, index) =>
        val blockInputs = {
          val inputSources: Vector[Map[String, CwlBlockInput]] = step.inputs.map { inp =>
            // sources of this step input - either a workflow input
            // like "file1" or a step output like "step1/file1"
            val sources = inp.source.map { src =>
              (src, IRParameter.encodeName(src.path.get))
            }
            val sourceParams = sources.foldLeft(Map.empty[String, Parameter]) {
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
              case (name, param) =>
                if (hasDefault || CwlOptional.isOptional(param.cwlType)) {
                  name -> OptionalBlockInput(name, param.cwlType)
                } else {
                  name -> RequiredBlockInput(name, param.cwlType)
                }
            }
          }
          inputSources.flatMap(_.toVector).groupBy(_._1).map {
            case (_, inputs) if inputs.toSet.size == 1 => inputs.head
            case (name, inputs) =>
              inputs.toSet
                .collectFirst {
                  case (_, inp: RequiredBlockInput) => name -> inp
                }
                .getOrElse(
                    throw new Exception(
                        s"multiple incompatible inputs ${inputs} map to name ${name}"
                    )
                )
          }
        }
        val procOutputs = step.run.outputs.collect {
          case o if o.hasName => o.name -> o
        }.toMap
        val blockOutputs = step.outputs.foldLeft(Map.empty[String, OutputParameter]) {
          case (accu, out) if procOutputs.contains(out.name) && !accu.contains(out.name) =>
            accu + (out.name -> procOutputs(out.name))
          case (_, out) =>
            throw new Exception(s"invalid or duplicate output parameter name ${out.name}")
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
