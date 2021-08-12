package dx.executor.cwl

import dx.api.{DxFile, InstanceTypeRequest}
import dx.core.Constants
import dx.cwl._
import dx.cwl.Document.Document
import dx.core.io.StreamFiles
import dx.core.ir.{Parameter, Type, Value}
import dx.core.languages.Language
import dx.core.languages.cwl.{CwlUtils, DxHintSchema, RequirementEvaluator, Target}
import dx.executor.{JobMeta, TaskExecutor}
import dx.util.{DockerUtils, FileUtils, JsUtils, TraceLevel}
import spray.json._

import java.net.URI
import java.nio.file.Files

object CwlTaskExecutor {
  def create(jobMeta: JobMeta,
             streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
             waitOnUpload: Boolean,
             checkInstanceType: Boolean): CwlTaskExecutor = {
    // when parsing a packed workflow as a String, we need to use a baseuri -
    // it doesn't matter what it is
    val parser = Parser.create(Some(URI.create("file:/null")), hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(jobMeta.sourceCode) match {
      case Some((version, "CommandLineTool" | "ExpressionTool" | "Workflow"))
          if Language.parse(version) == Language.CwlV1_2 =>
        ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL CommandLineTool, ExpressionTool,
               |or Workflow of a supported version:
               |${jobMeta.sourceCode}""".stripMargin
        )
    }
    val toolName = jobMeta.getExecutableAttribute("name") match {
      case Some(JsString(name)) => name
      case _                    => throw new Exception("missing executable name")
    }
    jobMeta.logger.trace(s"toolName: ${toolName}")
    // The main process may or may not have an ID. In the case it does not,
    // we supply a default ID fragment that is unlikely to be taken by the
    // user. Then, if the main process is a CommandLineTool or ExpressionTool,
    // we update the ID to use the tool name as the fragment; otherwise it is
    // a workflow and we look up the tool by name in the Map of processes nested
    // within the workflow.
    val defaultFrag = "___"
    val tool = parser.parseString(jobMeta.sourceCode, Some(defaultFrag)) match {
      case (tool: CommandLineTool, _) if tool.frag == defaultFrag =>
        tool.copy(id = Some(Identifier(namespace = None, frag = Some(toolName))))
      case (tool: CommandLineTool, _) => tool
      case (tool: ExpressionTool, _) if tool.frag == defaultFrag =>
        tool.copy(id = Some(Identifier(namespace = None, frag = Some(toolName))))
      case (tool: ExpressionTool, _) => tool
      case (_, doc: Document) =>
        doc.values.toVector.collect {
          case tool: CommandLineTool if tool.name == toolName => tool
          case tool: ExpressionTool if tool.name == toolName  => tool
        } match {
          case Vector(tool) => tool
          case Vector() =>
            throw new Exception(s"workflow does not contain a tool named ${toolName}")
          case v =>
            throw new Exception(s"more than one tool with name ${toolName}: ${v.mkString("\n")}")
        }
    }
    CwlTaskExecutor(tool, jobMeta, streamFiles, waitOnUpload, checkInstanceType)
  }
}

// TODO: add compile-time option to enable passing an overrides file to the
//  top-level workflow, and have the tool-specific overrides dispatched to
//  each task
// TODO: add compile-time option for --non-strict
// TODO: SHA1 checksums are computed for all outputs - we need to add these as
//  properties on the uploaded files so they can be propagated to downstream
//  CWL inputs
case class CwlTaskExecutor(tool: Process,
                           jobMeta: JobMeta,
                           streamFiles: StreamFiles.StreamFiles,
                           waitOnUpload: Boolean,
                           checkInstanceType: Boolean)
    extends TaskExecutor(jobMeta, streamFiles, waitOnUpload, checkInstanceType) {

  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger
  private val workerPaths = jobMeta.workerPaths

  override def executorName: String = "dxExecutorCwl"

  private lazy val inputParams: Map[String, InputParameter] = {
    tool.inputs.collect {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.flatMap(_.name).get -> param
    }.toMap
  }

  private lazy val runtime: Runtime = CwlUtils.createRuntime(workerPaths)

  private lazy val (cwlInputs: Map[String, (CwlType, CwlValue)], target: Option[String]) = {
    // CWL parameters can have '.' in their name
    val (irInputs, target) =
      jobMeta.primaryInputs.foldLeft(Map.empty[String, Value], Option.empty[String]) {
        case ((accu, None), (Target, Value.VString(targetName))) =>
          (accu, Some(targetName))
        case ((accu, target), (name, value)) =>
          (accu + (Parameter.decodeName(name) -> value), target)
      }
    val missingTypes = irInputs.keySet.diff(inputParams.keySet)
    if (missingTypes.nonEmpty) {
      throw new Exception(s"no type information given for input(s) ${missingTypes.mkString(",")}")
    }
    // convert IR to CWL values; discard auxiliary fields
    val evaluator = Evaluator.create(tool.requirements, tool.hints)
    val cwlInputs = inputParams.foldLeft(Map.empty[String, (CwlType, CwlValue)]) {
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
    if (logger.isVerbose) {
      val inputStr = tool.inputs
        .flatMap {
          case param if param.id.forall(_.name.forall(cwlInputs.contains)) =>
            val name = param.id.flatMap(_.name).get
            val (inputType, inputValue) = cwlInputs(name)
            Some(
                s"${name}: paramType=${param.cwlType}; inputType=${inputType}; inputValue=${inputValue}"
            )
          case other =>
            logger.trace(s"no input for parameter ${other}")
            None
        }
        .mkString("\n  ")
      logger.traceLimited(s"inputs:\n  ${inputStr}")
    }
    (cwlInputs, target)
  }

  private lazy val defaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = {
    CwlUtils.fromIRValues(jobMeta.defaultRuntimeAttrs, isInput = true)
  }

  private def getInstanceTypeRequest(
      inputs: Map[String, (CwlType, CwlValue)] = cwlInputs
  ): InstanceTypeRequest = {
    logger.traceLimited("calcInstanceType", minLevel = TraceLevel.VVerbose)
    val cwlEvaluator = Evaluator.create(tool.requirements, tool.hints)
    val ctx = CwlUtils.createEvaluatorContext(runtime)
    val env = cwlEvaluator.evaluateMap(inputs, ctx)
    val reqEvaluator = RequirementEvaluator(
        tool.requirements,
        tool.hints,
        env,
        workerPaths,
        tool.inputs.map(i => i.name -> i).toMap,
        defaultRuntimeAttrs
    )
    reqEvaluator.parseInstanceType
  }

  override protected lazy val getInstanceTypeRequest: InstanceTypeRequest = {
    getInstanceTypeRequest()
  }

  override protected def getInputsWithDefaults: Map[String, (Type, Value)] = {
    CwlUtils.toIR(cwlInputs)
  }

  override protected def streamFileForInput(parameterName: String): Boolean = {
    inputParams(parameterName).streamable
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
  ): (Map[String, (Type, Value)], Boolean, Option[Set[Int]], Set[Int]) = {
    val inputs = CwlUtils.fromIR(localizedInputs, typeAliases, isInput = true)
    val metaDir = workerPaths.getMetaDir(ensureExists = true)
    // write the CWL and input files
    val cwlPath = metaDir.resolve(s"tool.cwl")
    FileUtils.writeFileContent(cwlPath, jobMeta.sourceCode)
    val inputPath = metaDir.resolve(s"tool_input.json")
    val inputJson = CwlUtils.toJson(inputs)
    if (logger.isVerbose) {
      logger.trace(s"input JSON ${inputPath}:\n${inputJson.prettyPrint}")
    }
    JsUtils.jsToFile(inputJson, inputPath)
    // if a target is specified (a specific workflow step), add the
    // --single-process option
    val targetOpt = target.map(t => s"--single-process ${t}").getOrElse("")
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
    // TODO: remove --skip-schemas once we support them
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
         |    --skip-schemas \\
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
    val (successCodes, retryCodes) = tool match {
      case c: CommandLineTool =>
        (Some(c.successCodes), c.temporaryFailCodes)
      case _ => (Some(Set(0)), Set.empty[Int])
    }
    (localizedInputs, true, successCodes, retryCodes)
  }

  private lazy val outputParams: Map[String, OutputParameter] = {
    val outputParams = tool.outputs.map {
      case param if param.id.forall(_.name.isDefined) =>
        param.id.flatMap(_.name).get -> param
    }.toMap
    if (logger.isVerbose) {
      logger.traceLimited(s"outputParams=${outputParams}")
    }
    outputParams
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[String, (Type, Value)]
  ): (Map[String, (Type, Value)], Map[String, (Set[String], Map[String, String])]) = {
    // the outputs were written to stdout
    val stdoutFile = workerPaths.getStdoutFile()
    val localizedOutputs = if (Files.exists(stdoutFile)) {
      val allOutputs = JsUtils.jsFromFile(stdoutFile) match {
        case JsObject(outputs) => outputs
        case JsNull            => Map.empty[String, JsValue]
        case other             => throw new Exception(s"unexpected cwltool outputs ${other}")
      }
      val cwlOutputs = outputParams.map {
        case (name, param: WorkflowOutputParameter) if param.sources.nonEmpty =>
          val sources = param.sources.map(id => allOutputs.getOrElse(id.name.get, JsNull))
          val isArray = param.cwlType match {
            case _: CwlArray => true
            case CwlMulti(types) =>
              types.exists {
                case _: CwlArray => true
                case _           => false
              }
            case _ => false
          }
          val pickedAndMerged = if (sources.size == 1) {
            sources.head
          } else {
            val mergedValues = if (isArray) {
              param.linkMerge match {
                case LinkMergeMethod.MergeNested => sources
                case LinkMergeMethod.MergeFlattened =>
                  sources.flatMap {
                    case JsArray(items) => items
                    case value          => Vector(value)
                  }
              }
            } else {
              sources
            }
            val pickedValues: Vector[JsValue] = if (param.pickValue.nonEmpty) {
              val nonNull = mergedValues.filterNot(_ == JsNull)
              param.pickValue.get match {
                case PickValueMethod.FirstNonNull =>
                  Vector(
                      nonNull.headOption
                        .getOrElse(
                            throw new Exception(
                                s"all source values are null for parameter ${param}"
                            )
                        )
                  )
                case PickValueMethod.TheOnlyNonNull =>
                  if (nonNull.size == 1) {
                    Vector(nonNull.head)
                  } else {
                    throw new Exception(
                        s"there is not exactly one non-null value for parameter ${param}"
                    )
                  }
                case PickValueMethod.AllNonNull => nonNull
              }
            } else {
              mergedValues
            }
            if (isArray) {
              JsArray(pickedValues)
            } else if (pickedValues.size == 1) {
              pickedValues.head
            } else if (pickedValues.size > 1) {
              throw new Exception(
                  s"multiple output sources for non-array parameter ${param} that does not specify pickValue"
              )
            } else {
              JsNull
            }
          }
          name -> CwlValue.deserialize(pickedAndMerged, param.cwlType, typeAliases)
        case (name, param) if allOutputs.contains(name) =>
          name -> CwlValue.deserialize(allOutputs(name), param.cwlType, typeAliases)
        case (name, param) if CwlOptional.isOptional(param.cwlType) =>
          name -> (param.cwlType, NullValue)
        case (_, param) =>
          throw new Exception(s"missing value for output parameter ${param}")
      }
      CwlUtils.toIR(cwlOutputs)
    } else {
      Map.empty[String, (Type, Value)]
    }
    (localizedOutputs, Map.empty)
  }

  override protected def outputTypes: Map[String, Type] = {
    outputParams.map {
      case (name, param) => name -> CwlUtils.toIRType(param.cwlType)
    }
  }
}
