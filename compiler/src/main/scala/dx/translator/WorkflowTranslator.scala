package dx.translator

import dx.api.DxWorkflowStage
import dx.core.Constants
import dx.core.Constants.ReorgStatus
import dx.core.ir.RunSpec.{DefaultInstanceType, NoImage}
import dx.core.ir.Type.TFile
import dx.core.ir.Value.VFile
import dx.core.ir.{
  Application,
  Callable,
  ExecutableKindWfInputs,
  ExecutableKindWorkflowCustomReorg,
  ExecutableKindWorkflowOutputReorg,
  Parameter,
  SourceCode,
  Stage,
  StageInput,
  StaticInput,
  Workflow
}
import dx.util.Logger

abstract class WorkflowTranslator(wfName: String,
                                  availableDependencies: Map[String, Callable],
                                  reorgAttrs: ReorgSettings,
                                  logger: Logger) {
  type LinkedVar = (Parameter, StageInput)

  protected def isLocked: Boolean

  protected def standAloneWorkflow(setTarget: Boolean): SourceCode

  // generate a stage Id, this is a string of the form: 'stage-xxx'
  protected val fragNumIter: Iterator[Int] = Iterator.from(0)

  protected def getStageId(stageName: Option[String] = None): String = {
    stageName.map(name => s"stage-${name}").getOrElse(s"stage-${fragNumIter.next()}")
  }

  protected def getStage(stageName: Option[String] = None): DxWorkflowStage = {
    DxWorkflowStage(getStageId(stageName))
  }

  def translate: (Workflow, Vector[Callable], Vector[LinkedVar])

  /**
    * Create a preliminary applet to handle workflow input/outputs. This is
    * used only in the absence of workflow-level inputs/outputs.
    * @param wfName the workflow name
    * @param appletInputs the common applet inputs
    * @param stageInputs the common stage inputs
    * @param outputs the outputs
    * @return
    */
  protected def createCommonApplet(wfName: String,
                                   appletInputs: Vector[Parameter],
                                   stageInputs: Vector[StageInput],
                                   outputs: Vector[Parameter],
                                   blockPath: Vector[Int] = Vector.empty): (Stage, Application) = {
    val applet = Application(
        s"${wfName}_${Constants.CommonStage}",
        appletInputs,
        outputs,
        DefaultInstanceType,
        NoImage,
        ExecutableKindWfInputs(blockPath),
        standAloneWorkflow(setTarget = false)
    )
    logger.trace(s"Compiling common applet ${applet.name}")
    val stage = Stage(
        Constants.CommonStage,
        getStage(Some(Constants.CommonStage)),
        applet.name,
        stageInputs,
        outputs
    )
    (stage, applet)
  }

  /**
    * Creates an applet to reorganize the output files. We want to
    * move the intermediate results to a subdirectory.  The applet
    * needs to process all the workflow outputs, to find the files
    * that belong to the final results.
    * @param wfName workflow name
    * @param wfOutputs workflow outputs
    * @return
    */
  private def createReorgStage(wfName: String,
                               wfOutputs: Vector[LinkedVar]): (Stage, Application) = {
    val applet = Application(
        s"${wfName}_${Constants.ReorgStage}",
        wfOutputs.map(_._1),
        Vector.empty,
        DefaultInstanceType,
        NoImage,
        ExecutableKindWorkflowOutputReorg,
        standAloneWorkflow(setTarget = false)
    )
    logger.trace(s"Creating output reorganization applet ${applet.name}")
    // Link to the X.y original variables
    val inputs: Vector[StageInput] = wfOutputs.map(_._2)
    val stage =
      Stage(Constants.ReorgStage,
            getStage(Some(Constants.ReorgStage)),
            applet.name,
            inputs,
            Vector.empty[Parameter])
    (stage, applet)
  }

  private def createCustomReorgStage(wfOutputs: Vector[LinkedVar],
                                     appletId: String,
                                     reorgConfigFile: Option[String]): (Stage, Application) = {
    logger.trace(s"Creating custom output reorganization applet ${appletId}")
    val (statusParam, statusStageInput): LinkedVar = wfOutputs.filter {
      case (x, _) => x.name == ReorgStatus
    } match {
      case Vector(lvar) => lvar
      case other =>
        throw new Exception(
            s"Expected exactly one output with name ${ReorgStatus}, found ${other}"
        )
    }
    val configFile: Option[VFile] = reorgConfigFile.map(VFile(_))
    val appInputs = Vector(
        statusParam,
        Parameter(Constants.ReorgConfig, TFile, configFile)
    )
    val appletKind = ExecutableKindWorkflowCustomReorg(appletId)
    val applet = Application(
        appletId,
        appInputs,
        Vector.empty,
        DefaultInstanceType,
        NoImage,
        appletKind,
        standAloneWorkflow(setTarget = false)
    )
    // Link to the X.y original variables
    val inputs: Vector[StageInput] = configFile match {
      case Some(x) => Vector(statusStageInput, StaticInput(x))
      case _       => Vector(statusStageInput)
    }
    val stage =
      Stage(Constants.ReorgStage,
            getStage(Some(Constants.ReorgStage)),
            applet.name,
            inputs,
            Vector.empty[Parameter])
    (stage, applet)
  }

  def apply: Vector[Callable] = {
    val (irWf, irCallables, irOutputs) = translate
    // add a reorg applet if necessary
    val (updatedWf, updatedCallables) = reorgAttrs match {
      case DefaultReorgSettings(true) =>
        val (reorgStage, reorgApl) = createReorgStage(wfName, irOutputs)
        (irWf.copy(stages = irWf.stages :+ reorgStage), irCallables :+ reorgApl)
      case CustomReorgSettings(appUri, reorgConfigFile, true) if !isLocked =>
        val (reorgStage, reorgApl) = createCustomReorgStage(irOutputs, appUri, reorgConfigFile)
        (irWf.copy(stages = irWf.stages :+ reorgStage), irCallables :+ reorgApl)
      case _ =>
        (irWf, irCallables)
    }
    // validate workflow stages
    val allCallableNames = updatedCallables.map(_.name).toSet ++ availableDependencies.keySet
    val invalidStages =
      updatedWf.stages.filterNot(stage => allCallableNames.contains(stage.calleeName))
    if (invalidStages.nonEmpty) {
      val invalidStageDesc = invalidStages
        .map(stage => s"$stage.description, $stage.id.getId -> ${stage.calleeName}")
        .mkString("    ")
      throw new Exception(s"""|One or more stages reference missing tasks:
                              |stages: $invalidStageDesc
                              |callables: $allCallableNames
                              |""".stripMargin)
    }
    updatedCallables :+ updatedWf
  }
}
