package dx.executor.cwl

import dx.api.{DxFile, InstanceTypeRequest}
import dx.core.Constants
import dx.cwl._
import dx.core.io.StreamFiles.StreamFiles
import dx.core.ir.{Type, Value}
import dx.core.languages.cwl.{CwlUtils, RequirementEvaluator}
import dx.executor.{FileUploader, JobMeta, TaskExecutor}
import dx.util.{DockerUtils, FileUtils, JsUtils, TraceLevel}
import spray.json._

import java.nio.file.Files

object CwlTaskExecutor {
  def create(jobMeta: JobMeta,
             fileUploader: FileUploader,
             streamFiles: StreamFiles): CwlTaskExecutor = {
    if (!Parser.canParse(jobMeta.sourceCode)) {
      throw new Exception(
          s"""source code does not appear to be a CWL document of a supported version
             |${jobMeta.sourceCode}""".stripMargin
      )
    }
    val toolName = jobMeta.getExecutableAttribute("name") match {
      case Some(JsString(name)) => name
      case _                    => throw new Exception("missing executable name")
    }
    val tool =
      Parser.parseString(jobMeta.sourceCode, name = Some(toolName)) match {
        case tool: CommandLineTool => tool
        case other =>
          throw new Exception(s"expected CWL document to contain a CommandLineTool, not ${other}")
      }
    CwlTaskExecutor(tool, jobMeta, fileUploader, streamFiles)
  }
}

// TODO: add compile-time option to enable passing an overrides file to the
//  top-level workflow, and have the tool-specific overrides dispatched to
//  each task
// TODO: add compile-time option for --non-strict
// TODO: remove --skip-schemas once schema parsing is working
// TODO: SHA1 checksums are computed for all outputs - we need to add these as
//  properties on the uploaded files so they can be propagated to downstream
//  CWL inputs
case class CwlTaskExecutor(tool: CommandLineTool,
                           jobMeta: JobMeta,
                           fileUploader: FileUploader,
                           streamFiles: StreamFiles)
    extends TaskExecutor(jobMeta, fileUploader, streamFiles) {

  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger

  override def executorName: String = "dxExecutorCwl"

  private lazy val inputParams: Map[String, CommandInputParameter] = {
    tool.inputs.collect {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.get.name.get -> param
    }.toMap
  }

  private lazy val runtime: Runtime = CwlUtils.createRuntime(jobMeta.workerPaths)

  private def cwlInputs: Map[String, (CwlType, CwlValue)] = {
    val missingTypes = jobMeta.primaryInputs.keySet.diff(inputParams.keySet)
    if (missingTypes.nonEmpty) {
      throw new Exception(s"no type information given for input(s) ${missingTypes.mkString(",")}")
    }
    // convert IR to CWL values; discard auxiliary fields
    val evaluator = Evaluator.create(tool.requirements)
    inputParams.foldLeft(Map.empty[String, (CwlType, CwlValue)]) {
      case (env, (name, param)) =>
        val cwlTypes = param.types
        val (cwlType, cwlValue) = jobMeta.primaryInputs.get(name) match {
          case Some(irValue) =>
            CwlUtils.fromIRValue(irValue, cwlTypes, name)
          case None if param.default.isDefined =>
            val ctx = CwlUtils.createEvauatorContext(runtime, env)
            evaluator.evaluate(param.default.get, cwlTypes, ctx)
          case None if cwlTypes.exists(CwlOptional.isOptional) =>
            (CwlNull, NullValue)
          case _ =>
            throw new Exception(s"Missing required input ${name} to tool ${tool.id}")
        }
        env + (name -> (cwlType, cwlValue))
    }
  }

  private def printInputs(inputs: Map[String, (CwlType, CwlValue)]): Unit = {
    if (logger.isVerbose) {
      val inputStr = tool.inputs
        .flatMap {
          case param if param.id.forall(_.name.forall(inputs.contains)) =>
            val name = param.id.get.name.get
            Some(s"${name} -> (${param.types}, ${inputs.get(name)})")
          case other =>
            logger.trace(s"no input for parameter ${other}")
            None
        }
        .mkString("\n")
      logger.traceLimited(s"inputs: ${inputStr}")
    }
  }

  private lazy val defaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = {
    CwlUtils.fromIRValues(jobMeta.defaultRuntimeAttrs)
  }

  private def getRequiredInstanceTypeRequest(
      inputs: Map[String, (CwlType, CwlValue)] = cwlInputs
  ): InstanceTypeRequest = {
    logger.traceLimited("calcInstanceType", minLevel = TraceLevel.VVerbose)
    printInputs(inputs)
    val cwlEvaluator = Evaluator.create(tool.requirements)
    val ctx = CwlUtils.createEvauatorContext(runtime)
    val env = cwlEvaluator.evaluateMap(inputs, ctx)
    val reqEvaluator = RequirementEvaluator(
        tool.requirements,
        env,
        jobMeta.workerPaths,
        defaultRuntimeAttrs
    )
    reqEvaluator.parseInstanceType
  }

  override protected lazy val getRequiredInstanceTypeRequest: InstanceTypeRequest =
    getRequiredInstanceTypeRequest()

  override protected def getInputsWithDefaults: Map[String, (Type, Value)] = {
    val inputs = cwlInputs
    printInputs(inputs)
    CwlUtils.toIR(inputs)
  }

  override protected def streamFileForInput(parameterName: String): Boolean = {
    inputParams(parameterName).streamable.getOrElse(false)
  }

  private lazy val typeAliases: Map[String, CwlSchema] = {
    RequirementUtils.getSchemaDefs(tool.requirements)
  }

  override protected def getSchemas: Map[String, Type.TSchema] = {
    typeAliases.collect {
      case (name, record: CwlRecord) => name -> CwlUtils.toIRSchema(record)
    }
  }

  override protected def writeCommandScript(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    val inputs = CwlUtils.fromIR(localizedInputs, typeAliases)
    printInputs(inputs)
    val metaDir = jobMeta.workerPaths.getMetaDir(ensureExists = true)
    // write the CWL and input files
    val filePrefix = tool.id.name.getOrElse("tool")
    val cwlPath = metaDir.resolve(s"${filePrefix}.cwl")
    FileUtils.writeFileContent(cwlPath, jobMeta.sourceCode)
    val inputPath = metaDir.resolve(s"${filePrefix}_input.json")
    val inputJson = CwlUtils.toJson(inputs)
    JsUtils.jsToFile(inputJson, inputPath)
    // if a dx:// URI is specified for the Docker container, download it
    // and create an overrides file to override the value in the CWL file
    val overridesOpt = jobMeta
      .getExecutableDetail(Constants.DockerImage)
      .map { dockerImageJs =>
        val dockerImageDxFile = DxFile.fromJson(dxApi, dockerImageJs)
        val dockerUtils = DockerUtils(jobMeta.fileResolver, jobMeta.logger)
        val imageName = dockerUtils.getImage(dockerImageDxFile.asUri)
        val overridesJs = JsObject(
            "cwltool:overrides" -> JsObject(
                cwlPath.toString -> JsObject(
                    "requirements" -> JsObject(
                        "DockerRequirement" -> JsObject(
                            "dockerImageId" -> JsString(imageName)
                        )
                    )
                )
            )
        )
        val overridesPath = metaDir.resolve("overrides.json")
        JsUtils.jsToFile(overridesJs, overridesPath)
        s"--overrides ${overridesPath.toString}"
      }
      .getOrElse("")
    val command =
      s"""#!/bin/bash
         |(
         |  cwltool \\
         |    --basedir ${jobMeta.workerPaths.getRootDir(ensureExists = true)} \\
         |    --outdir ${jobMeta.workerPaths.getOutputFilesDir(ensureExists = true)} \\
         |    --tmpdir-prefix ${jobMeta.workerPaths.getTempDir(ensureExists = true)} \\
         |    --enable-pull \\
         |    --move-outputs \\
         |    --rm-container \\
         |    --rm-tmpdir \\
         |    --skip-schemas \\
         |    ${overridesOpt} ${cwlPath.toString} ${inputPath.toString}
         |) \\
         |> >( tee ${jobMeta.workerPaths.getStdoutFile(ensureParentExists = true)} ) \\
         |2> >( tee ${jobMeta.workerPaths.getStderrFile(ensureParentExists = true)} >&2 )
         |
         |echo $$? > ${jobMeta.workerPaths.getReturnCodeFile(ensureParentExists = true)}
         |
         |# make sure the files are on stable storage
         |# before leaving. This helps with stdin and stdout
         |# characters that may be in the fifo queues.
         |sync
         |""".stripMargin
    FileUtils.writeFileContent(
        jobMeta.workerPaths.getCommandFile(ensureParentExists = true),
        command,
        makeExecutable = true
    )
    localizedInputs
  }

  private lazy val outputParams: Map[String, CommandOutputParameter] = {
    tool.outputs.map {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.get.name.get -> param
    }.toMap
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    // the outputs were written to stdout
    val stdoutFile = jobMeta.workerPaths.getStdoutFile()
    if (Files.exists(stdoutFile)) {
      val cwlOutputs: Map[String, (CwlType, CwlValue)] = JsUtils.jsFromFile(stdoutFile) match {
        case JsObject(outputs) =>
          outputs.map {
            case (name, jsValue) =>
              val cwlTypes = outputParams(name).types
              name -> CwlValue.deserialize(jsValue, cwlTypes, typeAliases)
          }
        case JsNull => Map.empty
        case other  => throw new Exception(s"unexpected cwltool outputs ${other}")
      }
      CwlUtils.toIR(cwlOutputs)
    } else {
      Map.empty
    }
  }

  override protected def outputTypes: Map[String, Type] = {
    outputParams.map {
      case (name, param) if param.types.size == 1 =>
        name -> CwlUtils.toIRType(param.types.head)
      case (_, param) =>
        throw new Exception(s"Multi-type outputs are not supported ${param}")
    }
  }
}
