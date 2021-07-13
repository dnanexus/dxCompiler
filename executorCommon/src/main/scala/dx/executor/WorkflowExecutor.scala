package dx.executor

import dx.AppInternalException
import dx.api.{DxAnalysis, DxApp, DxApplet, DxExecution, DxFile, DxWorkflow, Field, FolderContents}
import dx.core.getVersion
import dx.core.ir.Type.TSchema
import dx.core.Constants
import dx.core.ir.{Block, ExecutableLink, Parameter, ParameterLink, Type, TypeSerde, Value}
import spray.json._
import dx.util.{Enum, LocalFileSource, TraceLevel}
import dx.util.CollectionUtils.IterableOnceExtensions

object WorkflowAction extends Enum {
  type WorkflowAction = Value
  val Inputs, Outputs, OutputReorg, CustomReorgOutputs, Run, Continue, Collect = Value
}

object WorkflowExecutor {
  val MaxNumFilesMoveLimit = 1000
  val IntermediateResultsFolder = "intermediate"
  val SeqNumber = "seqNumber"
  val JobNameLengthLimit = 50
  val ParentsKey = "parents___"
}

abstract class WorkflowExecutor[B <: Block[B]](jobMeta: JobMeta, separateOutputs: Boolean) {
  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger
  private val seqNumIter = Iterator.from(1)

  protected def nextSeqNum: Int = seqNumIter.next()

  val executorName: String

  protected def typeAliases: Map[String, TSchema]

  protected lazy val execLinkInfo: Map[String, ExecutableLink] =
    jobMeta.getExecutableDetail("execLinkInfo") match {
      case Some(JsObject(fields)) =>
        fields.map {
          case (key, link) =>
            key -> ExecutableLink.deserialize(link, typeAliases, jobMeta.dxApi)
        }
      case None =>
        Map.empty
      case other =>
        throw new Exception(s"Bad value ${other}")
    }

  protected def getExecutableLink(name: String): ExecutableLink = {
    execLinkInfo.getOrElse(
        name,
        throw new AppInternalException(s"Could not find linking information for ${name}")
    )
  }

  protected lazy val fqnDictTypes: Map[String, Type] =
    jobMeta.getExecutableDetail(Constants.WfFragmentInputTypes) match {
      case Some(jsv) => TypeSerde.deserializeSpec(jsv, typeAliases)
      case other     => throw new Exception(s"Bad value ${other}")
    }

  private lazy val jobInputs: Map[String, (Type, Value)] = {
    jobMeta.jsInputs.collect {
      case (name, jsValue) if !name.endsWith(ParameterLink.FlatFilesSuffix) =>
        val fqn = Parameter.decodeDots(name)
        val irType = fqnDictTypes.getOrElse(
            fqn,
            throw new Exception(s"Did not find variable ${fqn} (${name}) in the block environment")
        )
        val irValue = jobMeta.inputDeserializer.deserializeInputWithType(jsValue, irType, fqn)
        fqn -> (irType, irValue)
    }
  }

  protected def evaluateInputs(jobInputs: Map[String, (Type, Value)]): Map[String, (Type, Value)]

  private def evaluateInputs(): Map[String, ParameterLink] = {
    if (logger.isVerbose) {
      logger.traceLimited(s"dxCompiler version: ${getVersion}")
      logger.traceLimited(s"Environment: ${jobInputs}")
      logger.traceLimited("Artificial applet for workflow inputs")
    }
    val inputs = evaluateInputs(jobInputs)
    jobMeta.createOutputLinks(inputs)
  }

  protected def evaluateOutputs(jobInputs: Map[String, (Type, Value)],
                                addReorgStatus: Boolean): Map[String, (Type, Value)]

  private def evaluateOutputs(addReorgStatus: Boolean = false): Map[String, ParameterLink] = {
    if (logger.isVerbose) {
      logger.traceLimited(s"dxCompiler version: ${getVersion}")
      logger.traceLimited(s"Environment: ${jobInputs}")
      logger.traceLimited("Evaluating workflow outputs")
    }
    val outputs = evaluateOutputs(jobInputs, addReorgStatus)
    jobMeta.createOutputLinks(outputs)
  }

  /**
    * Evaluates any expressions in workflow outputs. The job output folder is set
    * to the job name, which is the concatenation of `name` and `nameDetail`.
    * @param executableLink the executable to launch
    * @param name the job name
    * @param inputs the job inputs
    * @param nameDetail: a suffix to add to the job name
    * @param instanceType the instance type to use for the new job
    * @param folder optional output folder
    * @param prefixOutputs whether to prefix output parameter names with the call
    *                      name when using manifest
    */
  protected def launchJob(executableLink: ExecutableLink,
                          name: String,
                          inputs: Map[String, (Type, Value)],
                          nameDetail: Option[String] = None,
                          instanceType: Option[String] = None,
                          folder: Option[String] = None,
                          prefixOutputs: Boolean = false): (DxExecution, String) = {
    val jobName: String = nameDetail.map(hint => s"${name} ${hint}").getOrElse(name)
    val outputFolder = if (separateOutputs) {
      folder.orElse(Some(jobName)).map {
        case f if !f.startsWith("/") => s"/${f}"
        case f                       => f
      }
    } else {
      None
    }

    val prefix = if (prefixOutputs) Some(name) else None
    val callInputsJs = JsObject(jobMeta.prepareSubjobInputs(inputs, executableLink, prefix))
    logger.traceLimited(s"""launchJob ${name} with arguments:
                           |${callInputsJs.prettyPrint}""".stripMargin,
                        minLevel = TraceLevel.VVerbose)

    // We may need to run a collect subjob. Add the the sequence number to each invocation, so the
    // collect subjob will be able to put the results back together in the correct order.
    val seqNum: Int = nextSeqNum

    // If this is a task that specifies the instance type at runtime, launch it in the requested instance.
    val dxExecution = executableLink.dxExec match {
      case app: DxApp =>
        app.newRun(
            jobName,
            callInputsJs,
            instanceType = instanceType,
            details = Some(JsObject(WorkflowExecutor.SeqNumber -> JsNumber(seqNum))),
            delayWorkspaceDestruction = jobMeta.delayWorkspaceDestruction,
            folder = outputFolder
        )
      case applet: DxApplet =>
        applet.newRun(
            jobName,
            callInputsJs,
            instanceType = instanceType,
            details = Some(JsObject(WorkflowExecutor.SeqNumber -> JsNumber(seqNum))),
            delayWorkspaceDestruction = jobMeta.delayWorkspaceDestruction,
            folder = outputFolder
        )
      case workflow: DxWorkflow =>
        workflow.newRun(
            jobName,
            callInputsJs,
            Some(JsObject(WorkflowExecutor.SeqNumber -> JsNumber(seqNum))),
            jobMeta.delayWorkspaceDestruction,
            folder = outputFolder
        )
      case other =>
        throw new Exception(s"Unsupported executable ${other}")
    }
    (dxExecution, jobName)
  }

  protected trait BlockContext {
    def block: Block[B]

    def launch(): Map[String, ParameterLink]

    def continue(): Map[String, ParameterLink]

    def collect(): Map[String, ParameterLink]

    def prettyFormat(): String

    // private variables may create files that need to be uploaded
    protected def uploadPrivateVariablePaths(
        env: Map[String, (Type, Value)]
    ): Map[String, ParameterLink] = {
      val filesToUpload = env
        .flatMap {
          case (name, (irType, irValue)) =>
            extractOutputFiles(name, irValue, irType, jobMeta.fileResolver)
        }
        .map {
          case local: LocalFileSource =>
            // if using manifests, we need to upload the files directly to the project
            val dest = if (jobMeta.useManifests) {
              Some(s"${jobMeta.manifestFolder}/${local.canonicalPath.getFileName.toString}")
            } else {
              None
            }
            local.address -> FileUpload(local.canonicalPath, dest)
            local.address -> FileUpload(local.canonicalPath, dest)
          case other =>
            throw new Exception(s"unexpected non-local file ${other}")
        }
        .toMap
      // upload the files, and map their local paths to their remote URIs
      val delocalizedPathToUri =
        fileUploader.upload(filesToUpload.values.toSet, wait = waitOnUpload).map {
          case (path, dxFile) => path -> dxFile.asUri
        }
    }
  }

  protected def evaluateBlockInputs(jobInputs: Map[String, (Type, Value)]): BlockContext

  private def evaluateFragInputs(): BlockContext = {
    if (logger.isVerbose) {
      logger.traceLimited(s"dxCompiler version: ${getVersion}")
      logger.traceLimited(s"link info=${execLinkInfo}")
      logger.traceLimited(s"Environment: ${jobInputs}")
    }
    val blockCtx = evaluateBlockInputs(jobInputs)
    logger.traceLimited(
        s"""|Block ${jobMeta.blockPath} to execute:
            |${blockCtx.prettyFormat()}
            |""".stripMargin
    )
    blockCtx
  }

  private def getAnalysisOutputFiles(analysis: DxAnalysis): Map[String, DxFile] = {
    // the analysis is updated asynchronously, so we need to check the dependsOn field
    // to determine if it is still waiting on the output of a previous stage before we
    // can proceed with the reorg
    val desc = Iterator
      .continually(
          analysis.describeNoCache(Set(Field.Input, Field.Output, Field.DependsOn))
      )
      .collectFirstDefined { a =>
        a.dependsOn match {
          case Some(Vector(jobId)) if jobId == jobMeta.jobId => Some(a)
          case None | Some(Vector())                         => Some(a)
          case _ =>
            Thread.sleep(3000)
            None
        }
      }
      .get
    val fileOutputs: Set[DxFile] = DxFile.findFiles(dxApi, desc.output.get).toSet
    val fileInputs: Set[DxFile] = DxFile.findFiles(dxApi, desc.input.get).toSet
    val analysisFiles: Vector[DxFile] = (fileOutputs -- fileInputs).toVector
    if (logger.isVerbose) {
      logger.traceLimited(s"analysis has ${fileOutputs.size} output files")
      logger.traceLimited(s"analysis has ${fileInputs.size} input files")
      logger.traceLimited(s"analysis has ${analysisFiles.size} real outputs")
    }
    val filteredAnalysisFiles: Map[String, DxFile] =
      if (analysisFiles.size > WorkflowExecutor.MaxNumFilesMoveLimit) {
        logger.traceLimited(
            s"WARNING: Large number of outputs (${analysisFiles.size}), not moving objects"
        )
        Map.empty
      } else {
        logger.traceLimited("Checking timestamps")
        // Retain only files that were created AFTER the analysis started
        val describedFiles = dxApi.describeFilesBulk(analysisFiles)
        val analysisCreated: java.util.Date = desc.getCreationDate
        describedFiles.collect {
          case dxFile if dxFile.describe().getCreationDate.compareTo(analysisCreated) >= 0 =>
            dxFile.id -> dxFile
        }.toMap
      }
    logger.traceLimited(s"analysis has ${filteredAnalysisFiles.size} verified output files")
    filteredAnalysisFiles
  }

  private def reorganizeOutputsDefault(): Map[String, ParameterLink] = {
    logger.traceLimited(s"dxCompiler version: ${getVersion}")
    val analysis = jobMeta.analysis.get
    val analysisFiles = getAnalysisOutputFiles(analysis)
    if (analysisFiles.isEmpty) {
      logger.trace(s"no output files to reorganize for analysis ${analysis.id}")
    } else {
      val outputIds: Set[String] = jobMeta.jsInputs
        .collect {
          case (name, jsValue) if !name.endsWith(ParameterLink.FlatFilesSuffix) =>
            val fqn = Parameter.decodeDots(name)
            if (!fqnDictTypes.contains(fqn)) {
              throw new Exception(
                  s"Did not find variable ${fqn} (${name}) in the block environment"
              )
            }
            DxFile.findFiles(dxApi, jsValue)
        }
        .flatten
        .map(_.id)
        .toSet
      val (exportFiles, intermediateFiles) = analysisFiles.partition {
        case (id, _) => outputIds.contains(id)
      }
      val exportNames: Vector[String] = exportFiles.values.map(_.describe().name).toVector
      logger.traceLimited(s"exportFiles=${exportNames}")
      if (intermediateFiles.isEmpty) {
        logger.trace("no intermediate files to reorganize for analysis ${analysis.id}")
      } else {
        val outputFolder = analysis.describe().folder
        logger.traceLimited(
            s"proj=${jobMeta.projectDesc.name} outFolder=${outputFolder}"
        )
        val intermediateFolder = s"${outputFolder}/${WorkflowExecutor.IntermediateResultsFolder}"
        val project = jobMeta.project
        val folderContents: FolderContents = project.listFolder(outputFolder)
        if (!folderContents.subFolders.contains(intermediateFolder)) {
          logger.traceLimited(s"Creating intermediate results sub-folder ${intermediateFolder}")
          project.newFolder(intermediateFolder, parents = true)
        } else {
          logger.traceLimited(
              s"Intermediate results sub-folder ${intermediateFolder} already exists"
          )
        }
        project.moveObjects(intermediateFiles.values.toVector, intermediateFolder)
      }
    }
    // reorg app has no outputs
    Map.empty
  }

  def apply(action: WorkflowAction.WorkflowAction): (Map[String, ParameterLink], String) = {
    try {
      val outputs: Map[String, ParameterLink] = action match {
        case WorkflowAction.Inputs =>
          evaluateInputs()
        case WorkflowAction.Run =>
          evaluateFragInputs().launch()
        case WorkflowAction.Continue =>
          evaluateFragInputs().continue()
        case WorkflowAction.Collect =>
          evaluateFragInputs().collect()
        case WorkflowAction.Outputs =>
          evaluateOutputs()
        case WorkflowAction.CustomReorgOutputs =>
          evaluateOutputs(addReorgStatus = true)
        case WorkflowAction.OutputReorg =>
          reorganizeOutputsDefault()
        case _ =>
          throw new Exception(s"Illegal workflow fragment operation ${action}")
      }
      jobMeta.writeOutputLinks(outputs)
      (outputs, s"success ${action}")
    } catch {
      case e: Throwable =>
        jobMeta.error(e)
        throw e
    }
  }
}
