package dx.executor.cwl

import dx.api.{DxFile, InstanceTypeRequest}
import dx.core.Constants
import dx.cwl._
import dx.core.io.StreamFiles.StreamFiles
import dx.core.ir.{Parameter, Type, Value}
import dx.core.languages.Language
import dx.core.languages.cwl.{CwlUtils, DxHintSchema, RequirementEvaluator}
import dx.executor.{FileUploader, JobMeta, TaskExecutor}
import dx.util.{DockerUtils, FileUtils, JsUtils, TraceLevel}
import spray.json._

import java.nio.file.Files

object CwlTaskExecutor {
  def create(jobMeta: JobMeta,
             fileUploader: FileUploader,
             streamFiles: StreamFiles): CwlTaskExecutor = {
    val parser = Parser.create(hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(jobMeta.sourceCode) match {
      case Some((version, "CommandLineTool")) if Language.parse(version) == Language.CwlV1_2 => ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL CommandLineTool document of a supported version
               |${jobMeta.sourceCode}""".stripMargin
        )
    }
    val toolName = jobMeta.getExecutableAttribute("name") match {
      case Some(JsString(name)) => name
      case _                    => throw new Exception("missing executable name")
    }
    val tool =
      parser.parseString(jobMeta.sourceCode, name = Some(toolName)) match {
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
  private val workerPaths = jobMeta.workerPaths

  override def executorName: String = "dxExecutorCwl"

  private lazy val inputParams: Map[String, CommandInputParameter] = {
    tool.inputs.collect {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.get.unqualifiedName.get -> param
    }.toMap
  }

  private lazy val runtime: Runtime = CwlUtils.createRuntime(workerPaths)

  private def cwlInputs: Map[String, (CwlType, CwlValue)] = {
    // CWL parameters can have '.' in their name
    val irInputs = jobMeta.primaryInputs.map {
      case (name, value) => Parameter.decodeName(name) -> value
    }
    val missingTypes = irInputs.keySet.diff(inputParams.keySet)
    if (missingTypes.nonEmpty) {
      throw new Exception(s"no type information given for input(s) ${missingTypes.mkString(",")}")
    }
    // convert IR to CWL values; discard auxiliary fields
    val evaluator = Evaluator.create(tool.requirements, tool.hints)
    inputParams.foldLeft(Map.empty[String, (CwlType, CwlValue)]) {
      case (env, (name, param)) =>
        val (cwlType: CwlType, cwlValue: CwlValue) = irInputs.get(name) match {
          case Some(irValue) =>
            CwlUtils.fromIRValue(irValue, param.cwlType, name, isInput = true)
          case None if param.default.isDefined =>
            val ctx = CwlUtils.createEvaluatorContext(runtime, env)
            evaluator.evaluate(param.default.get, param.cwlType, ctx)
          case None if CwlOptional.isOptional(param.cwlType) =>
            (param.cwlType, NullValue)
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
          case param if param.id.forall(_.unqualifiedName.forall(inputs.contains)) =>
            val name = param.id.get.unqualifiedName.get
            Some(s"${name} -> (${param.cwlType}, ${inputs.get(name)})")
          case other =>
            logger.trace(s"no input for parameter ${other}")
            None
        }
        .mkString("\n  ")
      logger.traceLimited(s"inputs:\n  ${inputStr}")
    }
  }

  private lazy val defaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = {
    CwlUtils.fromIRValues(jobMeta.defaultRuntimeAttrs, isInput = true)
  }

  private def getInstanceTypeRequest(
      inputs: Map[String, (CwlType, CwlValue)] = cwlInputs
  ): InstanceTypeRequest = {
    logger.traceLimited("calcInstanceType", minLevel = TraceLevel.VVerbose)
    printInputs(inputs)
    val cwlEvaluator = Evaluator.create(tool.requirements, tool.hints)
    val ctx = CwlUtils.createEvaluatorContext(runtime)
    val env = cwlEvaluator.evaluateMap(inputs, ctx)
    val reqEvaluator = RequirementEvaluator(
        tool.requirements,
        tool.hints,
        env,
        workerPaths,
        defaultRuntimeAttrs
    )
    reqEvaluator.parseInstanceType
  }

  override protected lazy val getInstanceTypeRequest: InstanceTypeRequest = {
    getInstanceTypeRequest()
  }

  override protected def getInputsWithDefaults: Map[String, (Type, Value)] = {
    val inputs = cwlInputs
    printInputs(inputs)
    CwlUtils.toIR(inputs)
  }

  override protected def streamFileForInput(parameterName: String): Boolean = {
    inputParams(parameterName).streamable.getOrElse(false)
  }

  private lazy val typeAliases: Map[String, CwlSchema] = {
    HintUtils.getSchemaDefs(tool.requirements)
  }

  override protected def getSchemas: Map[String, Type.TSchema] = {
    typeAliases.collect {
      case (name, record: CwlRecord) => name -> CwlUtils.toIRSchema(record)
    }
  }

  override protected def writeCommandScript(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    val inputs = CwlUtils.fromIR(localizedInputs, typeAliases, isInput = true)
    printInputs(inputs)
    val metaDir = workerPaths.getMetaDir(ensureExists = true)
    // write the CWL and input files
    val filePrefix = tool.getName.getOrElse("tool")
    val cwlPath = metaDir.resolve(s"${filePrefix}.cwl")
    FileUtils.writeFileContent(cwlPath, jobMeta.sourceCode)
    val inputPath = metaDir.resolve(s"${filePrefix}_input.json")
    val inputJson = CwlUtils.toJson(inputs)
    if (logger.isVerbose) {
      logger.trace(s"input JSON ${inputPath}:\n${inputJson.prettyPrint}")
    }
    JsUtils.jsToFile(inputJson, inputPath)
    // if a target is specified (a specific workflow step), add the --target option
    val targetOpt = jobMeta.targets.map { targets =>
      targets.map(t => s"--target ${t}").mkString(" ")
    }
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
        if (logger.isVerbose) {
          logger.trace(s"overrides JSON ${overridesPath}:\n${overridesJs.prettyPrint}")
        }
        JsUtils.jsToFile(overridesJs, overridesPath)
        s"--overrides ${overridesPath.toString}"
      }
      .getOrElse("")
    val command =
      s"""#!/bin/bash
         |(
         |  cwltool \\
         |    --basedir ${workerPaths.getRootDir(ensureExists = true)} \\
         |    --outdir ${workerPaths.getOutputFilesDir(ensureExists = true)} \\
         |    --tmpdir-prefix ${workerPaths.getTempDir(ensureExists = true)} \\
         |    --enable-pull \\
         |    --move-outputs \\
         |    --rm-container \\
         |    --rm-tmpdir \\
         |    ${targetOpt} ${overridesOpt} ${cwlPath.toString} ${inputPath.toString}
         |) \\
         |> >( tee ${workerPaths.getStdoutFile(ensureParentExists = true)} ) \\
         |2> >( tee ${workerPaths.getStderrFile(ensureParentExists = true)} >&2 )
         |
         |echo $$? > ${workerPaths.getReturnCodeFile(ensureParentExists = true)}
         |
         |# make sure the files are on stable storage
         |# before leaving. This helps with stdin and stdout
         |# characters that may be in the fifo queues.
         |sync
         |""".stripMargin
    FileUtils.writeFileContent(
        workerPaths.getCommandFile(ensureParentExists = true),
        command,
        makeExecutable = true
    )
    localizedInputs
  }

  private lazy val outputParams: Map[String, CommandOutputParameter] = {
    tool.outputs.map {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.get.unqualifiedName.get -> param
    }.toMap
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): Map[String, (Type, Value)] = {
    // the outputs were written to stdout
    val stdoutFile = workerPaths.getStdoutFile()
    if (Files.exists(stdoutFile)) {
      val cwlOutputs: Map[String, (CwlType, CwlValue)] = JsUtils.jsFromFile(stdoutFile) match {
        case JsObject(outputs) =>
          outputs.map {
            case (name, jsValue) =>
              val cwlTypes = outputParams(name).cwlType
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
      case (name, param) => name -> CwlUtils.toIRType(param.cwlType)
    }
  }
}
