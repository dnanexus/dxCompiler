package dx.core.languages.cwl

import dx.core.ir.Block
import dx.core.ir.BlockKind
import dx.cwl.{
  CommandLineTool,
  CwlOptional,
  ExpressionTool,
  Workflow,
  WorkflowInputParameter,
  WorkflowOutputParameter,
  WorkflowStep
}

case class CwlBlock(index: Int,
                    inputs: Map[String, (WorkflowInputParameter, Boolean)],
                    outputs: Map[String, WorkflowOutputParameter],
                    steps: Vector[WorkflowStep])
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

  override def getSubBlock(index: Int): CwlBlock = {
    (kind, target) match {
      case (BlockKind.ConditionalOneCall | BlockKind.ScatterOneCall, Some(step)) =>
        step.run match {
          case _: CommandLineTool if index == 0 =>
            CwlBlock(0, inputs, outputs, Vector(step))
          case wf: Workflow =>
            val innerBlocks = CwlBlock.createBlocks(wf)
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
      case (name, (param, optional)) =>
        val types = if (optional) {
          param.types.map(CwlOptional.ensureOptional)
        } else {
          param.types.map(CwlOptional.unwrapOptional)
        }
        s"${CwlUtils.prettyFormatTypes(types)} ${name}"
    }
    val outputStr = outputs.map {
      case (name, param) =>
        s"${CwlUtils.prettyFormatTypes(param.types)} ${name}"
    }
    val bodyStr = steps.map(prettyFormatStep).mkString("\n")
    s"""Block(${index}, ${kind})
       |  Inputs: ${inputStr}
       |  Outputs: ${outputStr}
       |  Body: ${bodyStr}""".stripMargin
  }
}

object CwlBlock {
  def createBlocks(wf: Workflow): Vector[CwlBlock] = {
    // ExpressionTools are bundled with the next tool or workflow step.
    val (parts, _) =
      wf.steps.foldLeft(Vector.empty[Vector[WorkflowStep]], Vector.empty[WorkflowStep]) {
        case ((parts, curPart), step) =>
          step.run match {
            case _: ExpressionTool =>
              (parts, curPart :+ step)
            case _ if curPart.isEmpty =>
              (parts :+ Vector(step), curPart)
            case _ =>
              (parts :+ curPart, Vector.empty)
          }
      }

    val wfInputs = wf.inputs.collect {
      case i if i.hasName => i.name -> i
    }.toMap
    val wfOutputs = wf.outputs.collect {
      case o if o.hasName => o.name -> o
    }.toMap

    parts.zipWithIndex.map {
      case (v, index) =>
        val blockInputs = v.foldLeft(Map.empty[String, (WorkflowInputParameter, Boolean)]) {
          case (accu, step) =>
            step.inputs.foldLeft(accu) {
              case (accu, inp) =>
                val missing = inp.source.toSet.diff(wfInputs.keySet)
                if (missing.nonEmpty) {
                  throw new Exception(s"undefined workflow input(s) ${missing.mkString(",")}")
                }
                val hasDefault = inp.default.isDefined
                inp.source.foldLeft(accu) {
                  case (accu, name) =>
                    val wfInput = wfInputs(name)
                    val optional = hasDefault || wfInput.types.forall {
                      case _: CwlOptional => true
                      case _              => false
                    }
                    accu.get(name) match {
                      case Some((_, true)) if !optional =>
                        accu + (name -> (wfInput, false))
                      case None =>
                        accu + (name -> (wfInputs(name), optional))
                      case _ => accu
                    }
                }
            }
        }
        val blockOutputs = v.foldLeft(Map.empty[String, WorkflowOutputParameter]) {
          case (accu, step) =>
            step.outputs.foldLeft(accu) {
              case (accu, out) if wfOutputs.contains(out.name) && !accu.contains(out.name) =>
                accu + (out.name -> wfOutputs(out.name))
              case (_, out) =>
                throw new Exception(s"invalid or duplicate output parameter name ${out.name}")
            }
        }
        CwlBlock(index, blockInputs, blockOutputs, v)
    }
  }
}
