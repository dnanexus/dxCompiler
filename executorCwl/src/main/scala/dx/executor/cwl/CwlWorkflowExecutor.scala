package dx.executor.cwl

import dx.core.ir.{Block, Type, TypeSerde, Value, ValueSerde}
import dx.core.languages.Language
import dx.core.languages.cwl.CwlUtils.toIRSchema
import dx.core.languages.cwl.{CwlBlock, CwlUtils, DxHintSchema}
import dx.cwl.{CwlOptional, CwlRecord, HintUtils, Parser, Workflow}
import dx.executor.{BlockContext, JobMeta, WorkflowExecutor}
import spray.json.JsString

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
      case (name, schema: CwlRecord) => name -> toIRSchema(schema)
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
                     s"  ${CwlUtils.prettyFormatTypes(inp.types)} ${inp.name}"
                   }
                   .mkString("\n")}""".stripMargin
          )
        }
        jobInputs
      case path =>
        val block: CwlBlock = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), path)
        val inputTypes = block.inputs.map {
          case (name, (param, optional)) =>
            if (optional) {
              name -> param.types.map(CwlOptional.ensureOptional)
            } else {
              name -> param.types
            }
        }
        if (logger.isVerbose) {
          logger.trace(
              s"""input parameters:
                 |${inputTypes
                   .map {
                     case (name, cwlType) =>
                       s"  ${CwlUtils.prettyFormatTypes(cwlType)} ${name}"
                   }
                   .mkString("\n")}""".stripMargin
          )
        }
        jobInputs.collect {
          case (name, (_, v)) if inputTypes.contains(name) =>
            // convert to CWL and back to IR to ensure the type is correct
            val (cwlType, cwlValue) =
              CwlUtils.fromIRValue(v, inputTypes(name), name, isInput = true)
            name -> CwlUtils.toIRValue(cwlValue, cwlType)
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
    val outputsParams = jobMeta.blockPath match {
      case Vector() => workflow.outputs
      case path =>
        val block = Block.getSubBlockAt(CwlBlock.createBlocks(workflow), path)

    }
  }

  override protected def evaluateBlockInputs(
      jobInputs: Map[String, (Type, Value)]
  ): BlockContext[CwlBlock] = ???
}
