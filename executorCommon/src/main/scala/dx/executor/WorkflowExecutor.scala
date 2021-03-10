package dx.executor

import dx.AppInternalException
import dx.api.{
  DxAnalysis,
  DxApp,
  DxApplet,
  DxExecution,
  DxFile,
  DxObject,
  DxPath,
  DxWorkflow,
  Field,
  FolderContents
}
import dx.core.getVersion
import dx.core.ir.Type.{TFile, TOptional, TSchema}
import dx.core.Constants
import dx.core.ir.Value.{VArray, VFile}
import dx.core.ir.{
  Block,
  BlockKind,
  ExecutableLink,
  Manifest,
  Parameter,
  ParameterLink,
  Type,
  TypeSerde,
  Value,
  ValueSerde
}
import dx.util.protocols.DxFileSource
import dx.util.CollectionUtils.IterableOnceExtensions
import spray.json._
import dx.util.{Enum, JsUtils, LocalFileSource, TraceLevel}

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

  // this method is exposed for unit testing
  def getComplexScatterName(items: Iterator[Option[String]],
                            maxLength: Int = JobNameLengthLimit): String = {
    // Create a name by concatenating the initial elements of the array.
    // Limit the total size of the name.
    val (_, strings, hasMore) =
      items.foldLeftWhile((-1, Vector.empty[String], false))(_._1 < maxLength) {
        case ((length, strings, _), Some(s)) =>
          val newLength = length + s.length + 1
          if (newLength > maxLength) {
            (newLength, strings, true)
          } else {
            (newLength, strings :+ s, false)
          }
        case ((length, strings, _), None) =>
          (length, strings, false)
      }
    val itemStr = strings.mkString(",")
    if (hasMore) {
      s"${itemStr},..."
    } else {
      itemStr
    }
  }
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

  private lazy val jobInputs: Map[String, (Type, Value)] = jobMeta.jsInputs.collect {
    case (name, jsValue) if !name.endsWith(ParameterLink.FlatFilesSuffix) =>
      val fqn = Parameter.decodeName(name)
      val irType = fqnDictTypes.getOrElse(
          fqn,
          throw new Exception(s"Did not find variable ${fqn} (${name}) in the block environment")
      )
      val irValue = jobMeta.inputDeserializer.deserializeInputWithType(jsValue, irType)
      fqn -> (irType, irValue)
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

    def env: Map[String, (Type, Value)]

    protected def prepareBlockOutputs(
        outputs: Map[String, ParameterLink]
    ): Map[String, ParameterLink] = {
      if (jobMeta.useManifests && outputs.contains(Constants.OutputManifest)) {
        outputs
      } else {
        val outputNames = block.outputNames
        logger.traceLimited(
            s"""|processOutputs
                |  env = ${env.keys}
                |  fragResults = ${outputs.keys}
                |  exportedVars = ${outputNames}
                |""".stripMargin,
            minLevel = TraceLevel.VVerbose
        )
        val filteredEnv = env.view.filterKeys(outputNames.contains).toMap
        val inputLinks = jobMeta.createOutputLinks(filteredEnv, validate = false)
        val outputLink = outputs.view.filterKeys(outputNames.contains).toMap
        inputLinks ++ outputLink
      }
    }

    protected def getFileName(path: String): String = {
      jobMeta.fileResolver.resolve(path) match {
        case local: LocalFileSource => local.canonicalPath.getFileName.toString
        case dx: DxFileSource       => dx.dxFile.getName
        case other                  => other.toString
      }
    }

    protected def truncate(scatterName: String): String = {
      if (scatterName.length > WorkflowExecutor.JobNameLengthLimit) {
        s"${scatterName.substring(0, WorkflowExecutor.JobNameLengthLimit - 3)}..."
      } else {
        scatterName
      }
    }

    /**
      * stick the IDs of all the parent jobs and all the child jobs to exclude
      * (i.e. the continue/collect jobs) into details - we'll use these in the
      * collect step
      * @param nextStart index at which to start the next scatter
      * @return
      */
    protected def createSubjobDetails(nextStart: Option[Int] = None,
                                      outputShape: Option[Vector[Int]] = None): JsValue = {
      val parents = jobMeta.getJobDetail(WorkflowExecutor.ParentsKey) match {
        case Some(JsArray(array)) => array.map(JsUtils.getString(_))
        case _                    => Vector.empty
      }
      // add the current job to the list of parents
      val allParents = parents :+ jobMeta.jobId
      val details = Vector(
          Some(WorkflowExecutor.ParentsKey -> JsArray(allParents.map(JsString(_)))),
          nextStart.map(i => Constants.ContinueStart -> JsNumber(i))
      ).flatten.toMap
      JsObject(details)
    }

    protected def prepareScatterResults(dxSubJob: DxExecution): Map[String, ParameterLink]

    /**
      * Lauch a job to continue a large scatter.
      * @param childJobs child jobs on which the continue job will depend
      * @param nextStart the index at which to continue the scatter
      * @return
      */
    protected def launchScatterContinue(
        childJobs: Vector[DxExecution],
        nextStart: Int
    ): Map[String, ParameterLink] = {
      assert(childJobs.nonEmpty)
      // Run a sub-job with the "continue" entry point.
      // We need to provide the exact same inputs.
      val dxSubJob: DxExecution = dxApi.runSubJob(
          "continue",
          Some(jobMeta.instanceTypeDb.defaultInstanceType.name),
          JsObject(jobMeta.rawJsInputs),
          childJobs,
          jobMeta.delayWorkspaceDestruction,
          Some(s"continue_scatter($nextStart)"),
          Some(createSubjobDetails(Some(nextStart)))
      )
      prepareScatterResults(dxSubJob)
    }

    protected def launchScatterCollect(
        childJobs: Vector[DxExecution]
    ): Map[String, ParameterLink] = {
      assert(childJobs.nonEmpty)
      // Run a sub-job with the "collect" entry point.
      // We need to provide the exact same inputs.
      val dxSubJob: DxExecution = dxApi.runSubJob(
          "collect",
          Some(jobMeta.instanceTypeDb.defaultInstanceType.name),
          JsObject(jobMeta.rawJsInputs),
          childJobs,
          jobMeta.delayWorkspaceDestruction,
          Some(s"collect_scatter"),
          Some(createSubjobDetails())
      )
      prepareScatterResults(dxSubJob)
    }

    protected def launchCall(blockIndex: Int): Map[String, ParameterLink]

    protected def launchConditional(): Map[String, ParameterLink]

    protected def launchScatter(): Map[String, ParameterLink]

    protected case class ChildExecution(execName: String,
                                        seqNum: Int,
                                        outputs: Map[String, JsValue],
                                        exec: DxExecution)

    private def parseOneResult(value: JsValue, excludeIds: Set[String]): Option[ChildExecution] = {
      val fields = value.asJsObject.fields
      val (exec, desc) = fields.get("id") match {
        case Some(JsString(id)) if excludeIds.contains(id) =>
          logger.trace(s"Ignoring result for job ${id}")
          return None
        case Some(JsString(id)) if id.startsWith("job-") =>
          val job = dxApi.job(id)
          val desc = fields("describe").asJsObject
          (job, desc)
        case Some(JsString(id)) if id.startsWith("analysis-") =>
          val analysis = dxApi.analysis(id)
          val desc = fields("describe").asJsObject
          (analysis, desc)
        case Some(other) =>
          throw new Exception(s"malformed id field ${other.prettyPrint}")
        case None =>
          throw new Exception(s"field id not found in ${value.prettyPrint}")
      }
      logger.trace(s"parsing desc ${desc} for ${exec}")
      val (execName, details, output) =
        desc.getFields("executableName", "details", "output") match {
          case Seq(JsString(execName), JsObject(details), JsObject(output)) =>
            (execName, details, output)
        }
      val seqNum = details.get(WorkflowExecutor.SeqNumber) match {
        case Some(JsNumber(i)) => i.toIntExact
        case other             => throw new Exception(s"Invalid seqNumber ${other}")
      }
      Some(ChildExecution(execName, seqNum, output, exec))
    }

    private def submitRequest(
        parentJobId: Option[String],
        cursor: JsValue,
        excludeIds: Set[String],
        limit: Option[Int]
    ): (Vector[ChildExecution], JsValue) = {
      val parentField: Map[String, JsValue] = parentJobId match {
        case None     => Map.empty
        case Some(id) => Map("parentJob" -> JsString(id))
      }
      val cursorField: Map[String, JsValue] = cursor match {
        case JsNull      => Map.empty
        case cursorValue => Map("starting" -> cursorValue)
      }
      val limitField: Map[String, JsValue] = limit match {
        case None    => Map.empty
        case Some(i) => Map("limit" -> JsNumber(i))
      }
      val describeField: Map[String, JsValue] = Map(
          "describe" -> JsObject(
              "fields" -> DxObject
                .requestFields(Set(Field.Output, Field.ExecutableName, Field.Details))
          )
      )
      val response = dxApi.findExecutions(parentField ++ cursorField ++ limitField ++ describeField)
      val results: Vector[ChildExecution] =
        response.fields.get("results") match {
          case Some(JsArray(results)) =>
            results.flatMap(res => parseOneResult(res, excludeIds))
          case Some(other) =>
            throw new Exception(s"malformed results field ${other.prettyPrint}")
          case None =>
            throw new Exception(s"missing results field ${response}")
        }
      (results, response.fields("next"))
    }

    private def findChildExecutions(parentJobId: Option[String],
                                    excludeIds: Set[String],
                                    limit: Option[Int] = None): Vector[ChildExecution] = {
      Iterator
        .unfold[Vector[ChildExecution], Option[JsValue]](Some(JsNull)) {
          case None => None
          case Some(cursor: JsValue) =>
            submitRequest(parentJobId, cursor, excludeIds, limit) match {
              case (Vector(), _)     => None
              case (results, JsNull) => Some(results, None)
              case (results, next)   => Some(results, Some(next))
            }
        }
        .toVector
        .flatten
        .sortWith(_.seqNum < _.seqNum)
    }

    protected def aggregateScatterJobOutputs: (Option[String], Vector[Map[String, JsValue]]) = {
      val childExecutions = jobMeta.getJobDetail(WorkflowExecutor.ParentsKey) match {
        case Some(JsArray(array)) =>
          val parentJobIds = array.map(JsUtils.getString(_))
          val excludeJobIds = parentJobIds.toSet + jobMeta.jobId
          parentJobIds.flatMap { parentJobId =>
            findChildExecutions(Some(parentJobId), excludeJobIds)
          }
        case _ =>
          val parentJob = jobMeta.parentJob match {
            case Some(job) => job
            case None =>
              throw new Exception(s"Can't get parent job for $jobMeta.jobDesc")
          }
          findChildExecutions(Some(parentJob.id), Set(jobMeta.jobId))
      }
      logger.trace(s"childExecs=${childExecutions}")

      if (jobMeta.useManifests) {
        // each job has an output manifest - we need to download them all
        childExecutions
          .foldLeft(Option.empty[String], Vector.empty[Map[String, JsValue]]) {
            case ((id, childOutputs), childExec) =>
              val manifestFile = childExec.outputs.get(Constants.OutputManifest) match {
                case Some(fileObj: JsObject) if DxFile.isLinkJson(fileObj) =>
                  Some(DxFile.fromJson(dxApi, fileObj))
                case Some(JsString(uri)) if uri.startsWith(DxPath.DxUriPrefix) =>
                  Some(dxApi.resolveFile(uri))
                case None =>
                  // maybe the applet doesn't have any outputs
                  None
                case other =>
                  throw new Exception(s"invalid manifest file value ${other}")
              }
              val (manifestId, manifestValues) = manifestFile
                .map { dxFile =>
                  val manifestJson = new String(dxApi.downloadBytes(dxFile)).parseJson
                  val manifest = Manifest.parse(manifestJson)
                  (manifest.id, manifest.jsValues)
                }
                .getOrElse((None, Map.empty[String, JsValue]))
              // all scatter jobs manifests should have the same ID
              val newId = (id, manifestId) match {
                case (None, None)                         => None
                case (None, Some(id))                     => Some(id)
                case (Some(id), None)                     => Some(id)
                case (Some(id1), Some(id2)) if id1 == id2 => Some(id1)
                case (Some(id1), Some(id2)) =>
                  throw new Exception(
                      s"scatter job output manifests had different IDs: ${id1} != ${id2}"
                  )
              }
              (newId, childOutputs :+ manifestValues)
          }
      } else {
        (None, childExecutions.map(_.outputs))
      }
    }

    protected def createScatterOutputArray(
        childOutputs: Vector[Map[String, JsValue]],
        name: String,
        irType: Type
    ): Value = {
      val nameEncoded = Parameter.encodeName(name)
      val arrayValue = childOutputs.flatMap { outputs =>
        (irType, outputs.get(nameEncoded)) match {
          case (_, Some(jsValue)) =>
            Some(jobMeta.inputDeserializer.deserializeInputWithType(jsValue, irType))
          case (TOptional(_), None) =>
            None
          case (_, None) =>
            // Required output that is missing
            throw new Exception(s"missing required field <${name}> in results")
        }
      }
      VArray(arrayValue)
    }

    protected def getScatterOutputs(
        childOutputs: Vector[Map[String, JsValue]]
    ): Map[String, (Type, Value)]

    private def collectScatter(): Map[String, ParameterLink] = {
      val (manifestId, childOutputs) = aggregateScatterJobOutputs
      val arrayValues: Map[String, (Type, Value)] = getScatterOutputs(childOutputs)
      if (arrayValues.isEmpty) {
        Map.empty
      } else if (jobMeta.useManifests) {
        if (manifestId.isEmpty) {
          throw new Exception("missing manifest Id")
        }
        // upload the merged manifest file
        val outputJson = arrayValues.map {
          case (name, (t, v)) => Parameter.encodeName(name) -> ValueSerde.serializeWithType(v, t)
        }
        val manifest = Manifest(outputJson, id = manifestId)
        val destination = s"${jobMeta.manifestFolder}/${jobMeta.jobId}_output.manifest.json"
        val manifestDxFile = dxApi.uploadString(manifest.toJson.prettyPrint, destination)
        val outputValues = Map(
            Constants.OutputManifest -> (TFile, VFile(manifestDxFile.asUri))
        )
        jobMeta.createOutputLinks(outputValues, validate = false)
      } else {
        jobMeta.createOutputLinks(arrayValues, validate = false)
      }
    }

    def launch(): Map[String, ParameterLink] = {
      val outputs: Map[String, ParameterLink] = block.kind match {
        case BlockKind.CallWithSubexpressions | BlockKind.CallFragment =>
          launchCall(block.index)
        case BlockKind.ConditionalOneCall | BlockKind.ConditionalComplex =>
          launchConditional()
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          assert(jobMeta.scatterStart == 0)
          launchScatter()
        case BlockKind.ExpressionsOnly => Map.empty
        case BlockKind.CallDirect =>
          throw new RuntimeException("unreachable state")
      }
      prepareBlockOutputs(outputs)
    }

    def continue(): Map[String, ParameterLink] = {
      val outputs: Map[String, ParameterLink] = block.kind match {
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          assert(jobMeta.scatterStart > 0)
          launchScatter()
        case _ =>
          throw new RuntimeException(s"cannot continue non-scatter block ${block}")
      }
      prepareBlockOutputs(outputs)
    }

    def collect(): Map[String, ParameterLink] = {
      val outputs: Map[String, ParameterLink] = block.kind match {
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          collectScatter()
        case _ =>
          throw new RuntimeException(s"cannot continue non-scatter block ${block}")
      }
      prepareBlockOutputs(outputs)
    }

    def prettyFormat(): String = {
      val envStr = if (env.isEmpty) {
        " <empty>"
      } else {
        "\n" + env
          .map {
            case (name, (t, v)) =>
              s"  ${name}: ${TypeSerde.toString(t)} ${ValueSerde.toString(v)}"
          }
          .mkString("\n")
      }
      s"""${block.prettyFormat}
         |Env:${envStr}
         |""".stripMargin
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
    val desc = analysis.describe(Set(Field.Input, Field.Output))
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
            val fqn = Parameter.decodeName(name)
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
