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
  FileUpload,
  FolderContents
}
import dx.core.{Constants, getVersion}
import dx.core.ir.{
  Block,
  BlockKind,
  DxName,
  ExecutableLink,
  Manifest,
  ParameterLink,
  Type,
  TypeSerde,
  Value,
  ValueSerde
}
import dx.core.ir.Type.{TFile, TSchema}
import dx.core.ir.Value.{TransformHandler, VArray, VFile, VString, WalkHandler}
import dx.util.{Enum, JsUtils, LocalFileSource, TraceLevel}
import dx.util.CollectionUtils.IterableOnceExtensions
import dx.util.protocols.{DxFileSource, DxFolderSource}
import spray.json._

import java.nio.file.Files

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
  private val dxNameFactory = jobMeta.dxNameFactory
  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger
  private val seqNumIter = Iterator.from(1)

  protected def nextSeqNum: Int = seqNumIter.next()

  val executorName: String

  protected def typeAliases: Map[String, TSchema]

  protected lazy val execLinkInfo: Map[String, ExecutableLink] =
    jobMeta.getExecutableDetail(Constants.ExecLinkInfo) match {
      case Some(JsObject(fields)) =>
        fields.map {
          case (key, link) =>
            key -> jobMeta.executableLinkDeserializer(link, typeAliases)
        }
      case None  => Map.empty
      case other => throw new Exception(s"invalid ${Constants.ExecLinkInfo} value: ${other}")
    }

  protected def getExecutableLink(name: String): ExecutableLink = {
    execLinkInfo.getOrElse(
        name,
        throw new AppInternalException(s"Could not find linking information for ${name}")
    )
  }

  protected lazy val dxNameToType: Map[DxName, Type] =
    jobMeta.getExecutableDetail(Constants.WfFragmentInputTypes) match {
      case Some(jsv) =>
        TypeSerde.deserializeSpec(jsv, typeAliases).map {
          case (encodedName, t) => dxNameFactory.fromDecodedName(encodedName) -> t
        }
      case other => throw new Exception(s"Bad value ${other}")
    }

  private lazy val jobInputs: Map[DxName, (Type, Value)] = {
    jobMeta.jsInputs.collect {
      case (dxName, jsValue) if !dxName.suffix.exists(_.endsWith(Constants.FlatFilesSuffix)) =>
        val irType = dxNameToType.getOrElse(
            dxName,
            throw new Exception(s"Did not find variable ${dxName} in the block environment")
        )
        val irValue =
          jobMeta.inputDeserializer.deserializeInputWithType(jsValue, irType, dxName.decoded)
        dxName -> (irType, irValue)
    }
  }

  protected def evaluateInputs(jobInputs: Map[DxName, (Type, Value)]): Map[DxName, (Type, Value)]

  private def evaluateInputs(): Map[DxName, ParameterLink] = {
    if (logger.isVerbose) {
      logger.traceLimited(s"dxCompiler version: ${getVersion}")
      logger.traceLimited(s"Environment: ${jobInputs}")
      logger.traceLimited("Artificial applet for workflow inputs")
    }
    val inputs = evaluateInputs(jobInputs)
    jobMeta.createOutputLinks(inputs)
  }

  protected def evaluateOutputs(jobInputs: Map[DxName, (Type, Value)],
                                addReorgStatus: Boolean): Map[DxName, (Type, Value)]

  private def evaluateOutputs(addReorgStatus: Boolean = false): Map[DxName, ParameterLink] = {
    if (logger.isVerbose) {
      logger.traceLimited(s"dxCompiler version: ${getVersion}")
      logger.traceLimited(s"Environment: ${jobInputs}")
      logger.traceLimited("Evaluating workflow outputs")
    }
    val outputs = evaluateOutputs(jobInputs, addReorgStatus)
    if (logger.isVerbose) {
      logger.traceLimited(s"Outputs: ${outputs}")
    }
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
                          inputs: Map[DxName, (Type, Value)],
                          nameDetail: Option[String] = None,
                          instanceType: Option[String] = None,
                          folder: Option[String] = None,
                          prefixOutputs: Boolean = false): (DxExecution, String) = {
    val jobName: String = nameDetail.map(hint => s"${name} ${hint}").getOrElse(name)
    val outputFolder =
      Option.when(separateOutputs)(DxFolderSource.ensureEndsWithSlash(folder.getOrElse(jobName)))

    val prefix = if (prefixOutputs) Some(name) else None
    val callInputsJs = JsObject(jobMeta.prepareSubjobInputs(inputs, executableLink, prefix).map {
      case (dxName, jsv) => dxName.encoded -> jsv
    })
    logger.traceLimited(s"""launchJob ${name} with arguments:
                           |${callInputsJs.prettyPrint}""".stripMargin,
                        minLevel = TraceLevel.VVerbose)

    // We may need to run a collect subjob. Add the the sequence number to each invocation, so the
    // collect subjob will be able to put the results back together in the correct order.
    val seqNum: Int = nextSeqNum
    val details = JsObject(WorkflowExecutor.SeqNumber -> JsNumber(seqNum))

    // If this is a task that specifies the instance type at runtime, launch it in the requested instance.
    val dxExecution = executableLink.dxExec match {
      case app: DxApp =>
        app.newRun(
            jobName,
            callInputsJs,
            instanceType = instanceType,
            details = Some(details),
            delayWorkspaceDestruction = jobMeta.delayWorkspaceDestruction,
            folder = outputFolder
        )
      case applet: DxApplet =>
        applet.newRun(
            jobName,
            callInputsJs,
            instanceType = instanceType,
            details = Some(details),
            delayWorkspaceDestruction = jobMeta.delayWorkspaceDestruction,
            folder = outputFolder
        )
      case workflow: DxWorkflow =>
        workflow.newRun(
            jobName,
            callInputsJs,
            Some(details),
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

    def env: Map[DxName, (Type, Value)]

    protected def prepareBlockOutputs(
        outputs: Map[DxName, ParameterLink]
    ): Map[DxName, ParameterLink] = {
      if (jobMeta.useManifests && outputs.contains(Constants.OutputManifest)) {
        outputs
      } else {
        val outputNames = block.outputNames
        logger.traceLimited(
            s"""|prepareBlockOutputs
                |  env = ${env.keys}
                |  fragResults = ${outputs.keys}
                |  exportedVars = ${outputNames}
                |""".stripMargin,
            minLevel = TraceLevel.Verbose
        )
        // Create output links from the values in the env, which came either
        // from job inputs or from evaluating private variables. Ignore any
        // variable in env that will be overridden by on in outputs.
        val filteredEnv = env.filter {
          case (k, _) if outputs.contains(k)           => false
          case (k, _) if block.outputNames.contains(k) => true
          case _                                       => false
        }
        val envLinks = jobMeta.createOutputLinks(filteredEnv, validate = false)
        // create output links from the exposed outputs
        val outputLink = outputs.view.filterKeys(outputNames.contains).toMap
        envLinks ++ outputLink
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
      * Stick the IDs of all the parent jobs and all the child jobs to exclude
      * (i.e. the continue/collect jobs) into details - we'll use these in the
      * collect step.
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
          nextStart.map(i => Constants.ContinueStart -> JsNumber(i)),
          outputShape.map(v => Constants.OutputShape -> JsArray(v.map(JsNumber(_))))
      ).flatten.toMap
      JsObject(details)
    }

    protected def prepareScatterResults(dxSubJob: DxExecution): Map[DxName, ParameterLink]

    /**
      * Lauch a job to continue a large scatter.
      * @param childJobs child jobs on which the continue job will depend
      * @param nextStart the index at which to continue the scatter
      * @return
      */
    protected def launchScatterContinue(
        childJobs: Vector[DxExecution],
        nextStart: Int,
        outputShape: Option[Vector[Int]] = None
    ): Map[DxName, ParameterLink] = {
      if (childJobs.nonEmpty) {
        // Run a sub-job with the "continue" entry point.
        // We need to provide the exact same inputs.
        val dxSubJob: DxExecution = dxApi.runSubJob(
            "continue",
            Some(jobMeta.instanceTypeDb.defaultInstanceType.name),
            JsObject(jobMeta.rawJsInputs.map {
              case (dxName, jsv) => dxName.encoded -> jsv
            }),
            childJobs,
            jobMeta.delayWorkspaceDestruction,
            Some(s"continue_scatter($nextStart)"),
            Some(createSubjobDetails(Some(nextStart), outputShape))
        )
        prepareScatterResults(dxSubJob)
      } else {
        Map.empty
      }
    }

    protected def launchScatterCollect(
        childJobs: Vector[DxExecution],
        outputShape: Option[Vector[Int]] = None
    ): Map[DxName, ParameterLink] = {
      if (childJobs.nonEmpty) {
        // Run a sub-job with the "collect" entry point.
        // We need to provide the exact same inputs.
        val dxSubJob: DxExecution = dxApi.runSubJob(
            "collect",
            Some(jobMeta.instanceTypeDb.defaultInstanceType.name),
            JsObject(jobMeta.rawJsInputs.map {
              case (dxName, jsv) => dxName.encoded -> jsv
            }),
            childJobs,
            jobMeta.delayWorkspaceDestruction,
            Some(s"collect_scatter"),
            Some(createSubjobDetails(outputShape = outputShape))
        )
        prepareScatterResults(dxSubJob)
      } else {
        Map.empty
      }
    }

    protected def launchCall(blockIndex: Int): Map[DxName, ParameterLink]

    protected def launchConditional(): Map[DxName, ParameterLink]

    protected def launchScatter(): Map[DxName, ParameterLink]

    protected case class ChildExecution(execName: String,
                                        seqNum: Int,
                                        outputs: Map[DxName, JsValue],
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
          case Seq(JsString(execName), JsObject(details), JsObject(jsOutput)) =>
            val output: Map[DxName, JsValue] = jsOutput.map {
              case (name, jsv) => dxNameFactory.fromEncodedName(name) -> jsv
            }
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

    protected def aggregateScatterJobOutputs
        : (Option[String], Option[String], Vector[Map[DxName, JsValue]]) = {
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
          .foldLeft(Option.empty[String], Option.empty[String], Vector.empty[Map[DxName, JsValue]]) {
            case ((id, execName, childOutputs), childExec) =>
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
                  val manifest = Manifest.parse(manifestJson, dxNameFactory)
                  (manifest.id, manifest.jsValues)
                }
                .getOrElse((None, Map.empty[DxName, JsValue]))
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
              val newExecName = (execName, childExec.execName) match {
                case (None, execName)                                       => Some(execName)
                case (Some(execName1), execName2) if execName1 == execName2 => Some(execName1)
                case (Some(execName1), execName2) =>
                  throw new Exception(
                      s"child execution names do not match: ${execName1} != ${execName2}"
                  )
              }
              (newId, newExecName, childOutputs :+ manifestValues)
          }
      } else {
        (None, None, childExecutions.map(_.outputs))
      }
    }

    protected def createScatterOutputArray(
        childOutputs: Vector[Map[DxName, JsValue]],
        name: DxName,
        irType: Type,
        execName: Option[String]
    ): Value = {
      val encodedNames = Vector(Some(name), execName.map(name.pushDecodedNamespace)).flatten
      VArray(childOutputs.flatMap { outputs =>
        encodedNames
          .collectFirst {
            case name if outputs.contains(name) =>
              jobMeta.inputDeserializer.deserializeInputWithType(outputs(name),
                                                                 irType,
                                                                 name.decoded)
          }
          .orElse {
            Option.when(!Type.isOptional(irType))(
                throw new Exception(s"missing required field <${name}> in results")
            )
          }
      })
    }

    protected def getScatterOutputs(
        childOutputs: Vector[Map[DxName, JsValue]],
        execName: Option[String]
    ): Map[DxName, (Type, Value)]

    private def collectScatter(): Map[DxName, ParameterLink] = {
      val (manifestId, execName, childOutputs) = aggregateScatterJobOutputs
      val arrayValues: Map[DxName, (Type, Value)] = getScatterOutputs(childOutputs, execName)
      if (arrayValues.isEmpty) {
        Map.empty
      } else if (jobMeta.useManifests) {
        if (manifestId.isEmpty) {
          throw new Exception("missing manifest Id")
        }
        // upload the merged manifest file
        val outputJson = arrayValues.map {
          case (dxName, (t, v)) => dxName -> ValueSerde.serializeWithType(v, t)
        }
        val manifest = Manifest(outputJson, id = manifestId)
        val destination =
          s"${jobMeta.manifestProjectAndFolder}/${jobMeta.jobId}_output.manifest.json"
        val manifestDxFile = dxApi.uploadString(manifest.toJson().prettyPrint, destination)
        val outputValues = Map(
            Constants.OutputManifest -> (TFile, VFile(manifestDxFile.asUri))
        )
        jobMeta.createOutputLinks(outputValues, validate = false)
      } else {
        jobMeta.createOutputLinks(arrayValues, validate = false)
      }
    }

    private object FileExtractor extends WalkHandler[Map[String, FileUpload]] {
      private val defaultDestParent = if (jobMeta.useManifests) {
        DxFolderSource.ensureEndsWithSlash(jobMeta.manifestProjectAndFolder)
      } else {
        "/"
      }

      def handleFile(uri: String,
                     optional: Boolean,
                     ctx: Map[String, FileUpload]): Map[String, FileUpload] = {
        jobMeta.fileResolver.resolve(uri) match {
          case local: LocalFileSource if optional && !Files.exists(local.canonicalPath) =>
            // ignore optional, non-existent files
            ctx
          case local: LocalFileSource if !Files.exists(local.canonicalPath) =>
            throw new Exception(
                s"required output file does not exist at ${local.canonicalPath}"
            )
          case local: LocalFileSource =>
            val destination = s"${defaultDestParent}${local.canonicalPath.getFileName.toString}"
            ctx + (local.address -> FileUpload(local.canonicalPath, Some(destination)))
          case _ => ctx
        }
      }

      override def apply(value: Value,
                         t: Option[Type],
                         optional: Boolean,
                         ctx: Map[String, FileUpload]): Option[Map[String, FileUpload]] = {
        (t, value) match {
          case (Some(TFile) | None, f: VFile) => Some(handleFile(f.uri, optional, ctx))
          case (Some(TFile), VString(uri))    => Some(handleFile(uri, optional, ctx))
          case _                              => None
        }
      }

      def extractFiles(value: Value, t: Type): Map[String, FileUpload] = {
        Value.walk(value, Some(t), Map.empty[String, FileUpload], FileExtractor)
      }
    }

    // private variables may create files that need to be uploaded
    protected def uploadPrivateVariablePaths(
        env: Map[DxName, (Type, Value)]
    ): Map[DxName, ParameterLink] = {
      val filesToUpload = env.values.flatMap {
        case (irType, irValue) => FileExtractor.extractFiles(irValue, irType)
      }.toMap
      // If there are any local files that were created as part of running this
      // fragment, upload them and replace the local URI with the remote one.
      // Otherwise just pass through env unchanged.
      val delocalizedOutputs = if (filesToUpload.nonEmpty) {
        // upload files
        val localPathToRemoteUri = jobMeta.uploadFiles(filesToUpload.values).map {
          case (path, dxFile) => path -> dxFile.asUri
        }
        // Replace the local paths in the output values with URIs. This requires
        // two look-ups: first to get the absolute Path associated with the file
        // value (which may be relative or absolute), and second to get the URI
        // associated with the Path. Returns an Optional[String] because optional
        // outputs may be null.
        object PathTranslator extends TransformHandler {
          private val localUriToPath = filesToUpload.map {
            case (uri, upload) => uri -> upload.source
          }
          override def apply(value: Value, t: Option[Type], optional: Boolean): Option[Value] = {
            val uri = (t, value) match {
              case (_, f: VFile)               => Some(f.uri)
              case (Some(TFile), VString(uri)) => Some(uri)
              case _                           => None
            }
            uri.flatMap(localUriToPath.get).flatMap(localPathToRemoteUri.get).map(VFile(_))
          }
        }
        env.map {
          case (name, (t, v)) => name -> (t, Value.transform(v, Some(t), PathTranslator))
        }
      } else {
        env
      }
      jobMeta.createOutputLinks(delocalizedOutputs)
    }

    def launch(): Map[DxName, ParameterLink] = {
      val outputs: Map[DxName, ParameterLink] = block.kind match {
        case BlockKind.CallWithSubexpressions | BlockKind.CallFragment =>
          launchCall(block.index)
        case BlockKind.ConditionalOneCall | BlockKind.ConditionalComplex =>
          launchConditional()
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          assert(jobMeta.scatterStart == 0)
          launchScatter()
        case BlockKind.ExpressionsOnly =>
          // upload files associated with any private variables that are exposed as block outputs
          uploadPrivateVariablePaths(env.filter {
            case (k, _) if block.inputNames.contains(k)  => false
            case (k, _) if block.outputNames.contains(k) => true
            case _                                       => false
          })
        case _ => throw new RuntimeException("unreachable state")
      }
      prepareBlockOutputs(outputs)
    }

    def continue(): Map[DxName, ParameterLink] = {
      val outputs: Map[DxName, ParameterLink] = block.kind match {
        case BlockKind.ScatterOneCall | BlockKind.ScatterComplex =>
          assert(jobMeta.scatterStart > 0)
          launchScatter()
        case _ =>
          throw new RuntimeException(s"cannot continue non-scatter block ${block}")
      }
      prepareBlockOutputs(outputs)
    }

    def collect(): Map[DxName, ParameterLink] = {
      val outputs: Map[DxName, ParameterLink] = block.kind match {
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
            case (dxName, (t, v)) =>
              s"  ${dxName}: ${TypeSerde.toString(t)} ${ValueSerde.toString(v)}"
          }
          .mkString("\n")
      }
      s"""${block.prettyFormat}
         |Env:${envStr}
         |""".stripMargin
    }
  }

  protected def evaluateBlockInputs(jobInputs: Map[DxName, (Type, Value)]): BlockContext

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
        val describedFiles = dxApi.describeFilesBulk(analysisFiles, searchWorkspaceFirst = true)
        val analysisCreated: java.util.Date = desc.getCreationDate
        describedFiles.collect {
          case dxFile if dxFile.describe().getCreationDate.compareTo(analysisCreated) >= 0 =>
            dxFile.id -> dxFile
        }.toMap
      }
    logger.traceLimited(s"analysis has ${filteredAnalysisFiles.size} verified output files")
    filteredAnalysisFiles
  }

  private def reorganizeOutputsDefault(): Map[DxName, ParameterLink] = {
    logger.traceLimited(s"dxCompiler version: ${getVersion}")
    val analysis = jobMeta.analysis.get
    val analysisFiles = getAnalysisOutputFiles(analysis)
    if (analysisFiles.isEmpty) {
      logger.trace(s"no output files to reorganize for analysis ${analysis.id}")
    } else {
      val outputIds: Set[String] = jobMeta.jsInputs
        .collect {
          case (dxName, jsValue) if !dxName.suffix.exists(_.endsWith(Constants.FlatFilesSuffix)) =>
            if (!dxNameToType.contains(dxName)) {
              throw new Exception(s"Did not find variable ${dxName} in the block environment")
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

  def apply(action: WorkflowAction.WorkflowAction): (Map[DxName, ParameterLink], String) = {
    val outputs: Map[DxName, ParameterLink] = action match {
      case WorkflowAction.Inputs             => evaluateInputs()
      case WorkflowAction.Run                => evaluateFragInputs().launch()
      case WorkflowAction.Continue           => evaluateFragInputs().continue()
      case WorkflowAction.Collect            => evaluateFragInputs().collect()
      case WorkflowAction.Outputs            => evaluateOutputs()
      case WorkflowAction.CustomReorgOutputs => evaluateOutputs(addReorgStatus = true)
      case WorkflowAction.OutputReorg        => reorganizeOutputsDefault()
      case _ =>
        throw new Exception(s"Illegal workflow fragment operation ${action}")
    }
    logger.traceLimited(s"outputLinks:\n  ${ParameterLink.prettyFormat(outputs)}")
    jobMeta.writeOutputLinks(outputs)
    (outputs, s"success ${action}")
  }
}
