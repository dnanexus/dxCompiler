package dx.executor.cwl

import dx.api.{DxFile, InstanceTypeRequest}
import dx.core.Constants
import dx.cwl._
import dx.cwl.Document.Document
import dx.core.io.StreamFiles
import dx.core.ir.{DxName, Type, Value}
import dx.core.languages.Language
import dx.core.languages.cwl.{CwlDxName, CwlUtils, DxHintSchema, RequirementEvaluator, Target}
import dx.executor.{JobMeta, TaskExecutor}
import dx.util.protocols.{DxFileSource, DxFolderSource}
import dx.util.{DockerUtils, FileUtils, JsUtils, LocalFileSource, TraceLevel}
import spray.json._

import java.net.URI
import java.nio.file.Files

object CwlTaskExecutor {
  def parseString(cwl: String): ParserResult = {
    // when parsing a packed workflow as a String, we need to use a baseuri - it doesn't matter what
    val parser = Parser.create(Some(URI.create("file:/null")), hintSchemas = Vector(DxHintSchema))
    parser.detectVersionAndClass(cwl) match {
      case (version, Some("CommandLineTool" | "ExpressionTool" | "Workflow"))
          if Language.parse(version) == Language.CwlV1_2 =>
        ()
      case _ =>
        throw new Exception(
            s"""source code does not appear to be a CWL CommandLineTool, ExpressionTool,
               |or Workflow of a supported version:
               |${cwl}""".stripMargin
        )
    }
    parser.parseString(cwl)
  }

  def create(jobMeta: JobMeta,
             streamFiles: StreamFiles.StreamFiles = StreamFiles.PerFile,
             checkInstanceType: Boolean): CwlTaskExecutor = {
    val toolName = jobMeta.getExecutableDetail(Constants.OriginalName) match {
      case Some(JsString(name)) => name
      case _                    => throw new Exception("missing executable name")
    }
    jobMeta.logger.trace(s"toolName: ${toolName}")
    val (tool, root) = parseString(jobMeta.sourceCode) match {
      case ParserResult(Some(tool: CommandLineTool), _, _, _) => (tool, None)
      case ParserResult(Some(expr: ExpressionTool), _, _, _)  => (expr, None)
      case ParserResult(Some(wf: Workflow), doc: Document, _, _) =>
        doc.values
          .foldLeft(Set.empty[Process]) {
            case (accu, tool: CommandLineTool) if tool.name == toolName =>
              accu + CwlUtils.simplifyProcess(tool)
            case (accu, expr: ExpressionTool) if expr.name == toolName =>
              accu + CwlUtils.simplifyProcess(expr)
            case (accu, _) => accu
          }
          .toVector match {
          case Vector(tool) => (tool, Some(wf))
          case Vector() =>
            throw new Exception(s"CWL document does not contain a tool named ${toolName}")
          case v =>
            throw new Exception(
                s"CWL document contains more than one tool with name ${toolName}: ${v.mkString("\n")}"
            )
        }
      case _ =>
        throw new Exception(
            s"CWL document cannot be parsed properly"
        )
    }
    CwlTaskExecutor(tool, root, jobMeta, streamFiles, checkInstanceType)
  }
}

// TODO: add compile-time option for --non-strict
// TODO: SHA1 checksums are computed for all outputs - we need to add these as properties on the
//  uploaded files so they can be propagated to downstream CWL inputs
case class CwlTaskExecutor(tool: Process,
                           root: Option[Workflow],
                           jobMeta: JobMeta,
                           streamFiles: StreamFiles.StreamFiles,
                           checkInstanceType: Boolean)
    extends TaskExecutor(jobMeta, streamFiles, checkInstanceType) {

  private val dxApi = jobMeta.dxApi
  private val logger = jobMeta.logger
  private val workerPaths = jobMeta.workerPaths

  override def executorName: String = "dxExecutorCwl"

  private lazy val inputParams: Map[DxName, InputParameter] = {
    tool.inputs.collect {
      case param if param.id.isDefined =>
        CwlDxName.fromSourceName(param.name) -> param
    }.toMap
  }

  private lazy val runtime: Runtime = CwlUtils.createRuntime(workerPaths)

  private lazy val overridesJs: Map[String, JsValue] = jobMeta.jsOverrides match {
    case Some(JsObject(fields)) if fields.contains("requirements") || fields.contains("hints") =>
      fields
    case Some(requirements: JsArray) => Map("requirements" -> requirements)
    case None                        => Map.empty
    case other =>
      throw new Exception(s"invalid overrides ${other}")
  }

  // the requirements and hints for the tool, which may be augmented by runtime overrides
  private lazy val (requirements, hints) = {
    if (overridesJs.nonEmpty) {
      // create a minimal document and then parse it to get the AST values
      val doc = overridesJs ++ Map(
          "cwlVersion" -> JsString("v1.2"),
          "class" -> JsString("CommandLineTool"),
          "id" -> JsString("overrides"),
          "inputs" -> JsArray(),
          "outputs" -> JsArray()
      )
      CwlTaskExecutor.parseString(JsObject(doc).prettyPrint) match {
        case ParserResult(Some(overridesTool: CommandLineTool), _, _, _) =>
          (tool.requirements ++ overridesTool.requirements, tool.hints ++ overridesTool.hints)
        case other =>
          throw new Exception(s"expected CommandLineTool, not ${other}")
      }
    } else {
      (tool.requirements, tool.hints)
    }
  }

  private val targetRegex = "(?:(.+)#)?(.+)".r

  private lazy val (cwlInputs: Map[DxName, (CwlType, CwlValue)],
                    targetProcess: Option[String],
                    targetStep: Option[String]) = {
    val (irInputs, targetProcess, targetStep) =
      jobMeta.primaryInputs.foldLeft(Map.empty[DxName, Value],
                                     Option.empty[String],
                                     Option.empty[String]) {
        case ((accu, None, None), (Target, Value.VString(targetRegex(targetWf, targetStep)))) =>
          (accu, Option(targetWf), Some(targetStep))
        case ((accu, targetProcess, targetStep), (dxName, value)) =>
          (accu + (dxName -> value), targetProcess, targetStep)
      }
    val missingTypes = irInputs.keySet.diff(inputParams.keySet)
    if (missingTypes.nonEmpty) {
      throw new Exception(s"no type information given for input(s) ${missingTypes.mkString(",")}")
    }
    // if this is a step in a workflow, get default values from step inputs
    val stepDefaults = Option
      .when(targetStep.isDefined && root.isDefined) {
        root.get.steps.collectFirst {
          case step if step.id.map(_.frag).contains(targetStep.get) =>
            val stepDefaults: Map[DxName, CwlValue] = step.inputs.collect {
              case inp if inp.default.isDefined =>
                CwlDxName.fromSourceName(inp.name) -> inp.default.get
            }.toMap
            stepDefaults
        }
      }
      .flatten
      .getOrElse(Map.empty)
    // convert IR to CWL values; discard auxiliary fields
    val evaluator = Evaluator.create(requirements, hints)
    val cwlInputs = inputParams.foldLeft(Map.empty[DxName, (CwlType, CwlValue)]) {
      case (env, (name, param)) =>
        val (cwlType: CwlType, cwlValue: CwlValue) = irInputs.get(name) match {
          case Some(irValue) =>
            CwlUtils.fromIRValueWithType(irValue,
                                         param.cwlType,
                                         name.decoded,
                                         isInput = true,
                                         fileResolver = jobMeta.fileResolver)
          case None if stepDefaults.contains(name) =>
            val ctx = CwlUtils.createEvaluatorContext(env = env.map {
              case (dxName, tv) => dxName.decoded -> tv
            }, runtime = runtime)
            evaluator.evaluate(stepDefaults(name), param.cwlType, ctx, coerce = true)
          case None if param.default.isDefined =>
            val ctx = CwlUtils.createEvaluatorContext(env = env.map {
              case (dxName, tv) => dxName.decoded -> tv
            }, runtime = runtime)
            evaluator.evaluate(param.default.get, param.cwlType, ctx, coerce = true)
          case None if CwlOptional.isOptional(param.cwlType) =>
            (param.cwlType, NullValue)
          case _ =>
            throw new Exception(s"Missing required input ${name} to tool ${tool.name}")
        }
        env + (name -> (cwlType, cwlValue))
    }
    if (logger.isVerbose) {
      val inputStr = tool.inputs
        .flatMap { param =>
          param.id.map(_.name) match {
            case Some(name) =>
              val dxName = CwlDxName.fromSourceName(name)
              cwlInputs.get(dxName) match {
                case Some((inputType, inputValue)) =>
                  Some(
                      s"${name}: paramType=${param.cwlType}; inputType=${inputType}; inputValue=${inputValue}"
                  )
                case _ =>
                  logger.trace(s"no input for parameter ${param}")
                  None
              }
            case _ =>
              logger.trace(s"no name for parameter ${param}")
              None
          }
        }
        .mkString("\n  ")
      logger.traceLimited(s"inputs:\n  ${inputStr}")
    }
    (cwlInputs, targetProcess, targetStep)
  }

  override protected def getInputVariables: Map[DxName, (Type, Value)] = {
    CwlUtils.toIR(cwlInputs)
  }

  private lazy val typeAliases: Map[String, CwlSchema] = {
    HintUtils.getSchemaDefs(requirements)
  }

  private lazy val defaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = {
    CwlUtils.fromIRValues(jobMeta.defaultRuntimeAttrs,
                          isInput = true,
                          fileResolver = jobMeta.fileResolver)
  }

  override protected def getInstanceTypeRequest(
      inputs: Map[DxName, (Type, Value)]
  ): InstanceTypeRequest = {
    logger.traceLimited("calcInstanceType", minLevel = TraceLevel.VVerbose)
    val cwlEvaluator = Evaluator.create(requirements, hints)
    val cwlInputs =
      CwlUtils.fromIR(inputs, typeAliases, isInput = true, fileResolver = jobMeta.fileResolver)
    val ctx = CwlUtils.createEvaluatorContext(runtime = runtime)
    val env = cwlEvaluator.evaluateMap(cwlInputs.map {
      case (dxName, tv) => dxName.decoded -> tv
    }, ctx)
    val reqEvaluator = RequirementEvaluator(
        requirements,
        hints,
        env,
        workerPaths,
        tool.inputs.map(i => i.name -> i).toMap,
        defaultRuntimeAttrs,
        fileResolver = jobMeta.fileResolver,
        dxApi = jobMeta.dxApi
    )
    reqEvaluator.parseInstanceType
  }

  override protected def streamFileForInput(parameterName: DxName): Boolean = {
    inputParams(parameterName).streamable
  }

  override protected def getSchemas: Map[String, Type.TSchema] = {
    typeAliases.collect {
      case (name, record: CwlRecord) => name -> CwlUtils.toIRSchema(record)
    }
  }

  // scan through the source file and update any hard-coded files/directories
  // with a location in `dependencies`
  private def updateSourceCode(dependencies: Map[String, (Type, Value)],
                               finalOverrides: Map[String, JsValue],
                               targetToolID: String): String = {
    def updateListing(jsv: JsValue): JsValue = {
      jsv match {
        case JsArray(entries) =>
          JsArray(entries.map {
            case file: JsObject if file.fields.get("class").contains(JsString("File")) =>
              file.fields
                .get("location")
                .orElse(file.fields.get("path"))
                .flatMap {
                  case JsString(uri) =>
                    jobMeta.fileResolver.resolve(uri) match {
                      case DxFileSource(dxFile, _) =>
                        dependencies.get(dxFile.asUri).map {
                          case (Type.TFile, f: Value.VFile) =>
                            val (_, cwlValue) =
                              CwlUtils.fromIRValueWithType(f,
                                                           CwlFile,
                                                           "",
                                                           isInput = true,
                                                           fileResolver = jobMeta.fileResolver)
                            cwlValue.toJson
                          case other =>
                            throw new Exception(s"expected file value for uri ${uri}, not ${other}")
                        }
                      case _ => None
                    }
                  case other =>
                    throw new Exception(s"expected file location/path to be a string, not ${other}")
                }
                .getOrElse(file)
            case dir: JsObject if dir.fields.get("class").contains(JsString("Directory")) =>
              dir.fields
                .get("location")
                .orElse(dir.fields.get("path"))
                .flatMap {
                  case JsString(uri) =>
                    jobMeta.fileResolver.resolveDirectory(uri) match {
                      case fs: DxFolderSource => dependencies.get(fs.address)
                      case _                  => None
                    }
                  case other =>
                    throw new Exception(
                        s"expected directory location/path to be a string, not ${other}"
                    )
                } match {
                case Some((Type.TDirectory, f: Value.VFolder)) =>
                  val (_, cwlValue) =
                    CwlUtils.fromIRValueWithType(f,
                                                 CwlDirectory,
                                                 "",
                                                 isInput = true,
                                                 fileResolver = jobMeta.fileResolver)
                  cwlValue.toJson
                case Some(other) =>
                  throw new Exception(s"expected directory URI to be a VFolder, not ${other}")
                case None if dir.fields.contains("listing") =>
                  JsObject(dir.fields + ("listing" -> updateListing(dir.fields("listing"))))
                case None => dir
              }
            case arr: JsArray => updateListing(arr)
            case other        => other
          })
        case _ => jsv
      }
    }
    def inner(jsv: JsValue): JsValue = {
      jsv match {
        case JsObject(fields) =>
          JsObject(fields.map {
            case ("requirements", JsArray(reqs)) =>
              if (logger.isVerbose) {
                logger.trace(s"field id:${fields.get("id")}")
                logger.trace(s"req before override:\n${reqs}")
                logger.trace(s"targetToolID: ${targetToolID}")
                logger.trace(s"override reqs: ${JsUtils
                  .getOptionalValues(finalOverrides, "requirements")
                  .getOrElse(Vector.empty)}")
              }
              val finalreqs = Option
                .when(fields.get("id").equals(JsString(targetToolID)))(
                    reqs ++ JsUtils
                      .getOptionalValues(finalOverrides, "requirements")
                      .getOrElse(Vector.empty)
                )
                .getOrElse(reqs)
              if (logger.isVerbose) {
                logger.trace(s"req after override:\n${reqs}")
              }
              "requirements" -> JsArray(finalreqs.map {
                case JsObject(reqFields)
                    if reqFields.get("class").contains(JsString("InitialWorkDirRequirement")) =>
                  JsObject(reqFields + ("listing" -> updateListing(reqFields("listing"))))
                case other => other
              })
            case ("$schemas", JsArray(schemas)) =>
              "$schemas" -> JsArray(
                  schemas.map {
                    case JsString(uri) if dependencies.contains(uri) =>
                      val (_, v) = dependencies(uri)
                      v match {
                        case f: Value.VFile =>
                          jobMeta.fileResolver.resolve(f.uri) match {
                            case fs: LocalFileSource => JsString(fs.canonicalPath.toString)
                            case fs                  => JsString(fs.address)
                          }
                        case other =>
                          throw new Exception(
                              s"expected schema URI ${uri} to resolve to a file value, not ${other}"
                          )
                      }
                    case other => other
                  }
              )
            case (k, v) => k -> inner(v)
          })
        case JsArray(items) => JsArray(items.map(inner))
        case _              => jsv
      }
    }
    if (logger.isVerbose) {
      logger.trace(s"before inner")
    }
    val updatedSource = inner(jobMeta.sourceCode.parseJson).prettyPrint
    if (logger.isVerbose) {
      logger.trace(s"updated source code:\n${updatedSource}")
    }
    updatedSource
  }

  private lazy val rawStep: Option[String] = {
    targetStep.map(_.split('/').grouped(2).flatMap(_.take(1)).toArray.mkString("/"))
  }
  override protected def writeCommandScript(
      localizedInputs: Map[DxName, (Type, Value)],
      localizedDependencies: Option[Map[String, (Type, Value)]]
  ): (Boolean, Option[Set[Int]]) = {
    val inputs = CwlUtils.fromIR(localizedInputs,
                                 typeAliases,
                                 isInput = true,
                                 fileResolver = jobMeta.fileResolver)
    val metaDir = workerPaths.getMetaDir(ensureExists = true).asJavaPath
    // write the CWL and input files
    val cwlPath = metaDir.resolve(s"tool.cwl")
    val cwlPathStr = targetProcess.map(p => s"${cwlPath}#${p}").getOrElse(cwlPath.toString)
    val inputPath = metaDir.resolve("tool_input.json")
    val inputJson = CwlUtils.toJson(inputs)
    if (logger.isVerbose) {
      logger.trace(s"input JSON ${inputPath}:\n${inputJson.prettyPrint}")
    }
    JsUtils.jsToFile(inputJson, inputPath)
    // if a target is specified (a specific workflow step), add the --single-process option
    val targetOpt = rawStep.map(t => s"--single-process ${t}").getOrElse("")
    // if a dx:// URI is specified for the Docker container, download it and create an overrides
    // file to override the value in the CWL file
    val requirementOverrides =
      JsUtils.getOptionalValues(overridesJs, "requirements").getOrElse(Vector.empty)
    val finalOverrides = Option
      .when(
          !requirementOverrides.exists {
            case JsObject(fields) if fields.get("class").contains(JsString("DockerRequirement")) =>
              true
            case _ => false
          }
      )(
          jobMeta
            .getExecutableDetail(Constants.DockerImage)
            .map { dockerImageJs =>
              val dockerImageDxFile = DxFile.fromJson(dxApi, dockerImageJs)
              val dockerUtils = DockerUtils(jobMeta.fileResolver, jobMeta.logger)
              val imageName = dockerUtils.getImage(dockerImageDxFile.asUri)
              JsObject(
                  "class" -> JsString("DockerRequirement"),
                  "dockerImageId" -> JsString(imageName)
              )
            }
      )
      .flatten
      .map { dockerOverride =>
        overridesJs ++ Map(
            "requirements" -> JsArray(requirementOverrides ++ Vector(dockerOverride))
        )
      }
      .getOrElse(overridesJs)
    val overridesOpt = if (finalOverrides.nonEmpty) {
      val overridesDocJs = JsObject(
          "cwltool:overrides" -> JsObject(
              "#main" -> JsObject(finalOverrides),
              s"#${tool.frag}" -> JsObject(finalOverrides)
          )
      )
      val overridesDocPath = metaDir.resolve("overrides.json")
      if (logger.isVerbose) {
        logger.trace(s"overrides JSON ${overridesDocPath}:\n${overridesDocJs.prettyPrint}")
      }
      JsUtils.jsToFile(overridesDocJs, overridesDocPath)
      s"--overrides ${overridesDocPath.toString}"
    } else {
      ""
    }
    // update the source code if necessary
    val sourceCode = localizedDependencies
      .map(dep => updateSourceCode(dep, finalOverrides, tool.name))
      .getOrElse(jobMeta.sourceCode)
    if (logger.isVerbose) {
      logger.trace(
          s"""Executing CWL:
             |${sourceCode}""".stripMargin
      )
    }
    FileUtils.writeFileContent(cwlPath, sourceCode)
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
         |    ${targetOpt} ${overridesOpt} ${cwlPathStr} ${inputPath.toString}
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
        workerPaths.getCommandFile(ensureParentExists = true).asJavaPath,
        command,
        makeExecutable = true
    )
    // We are testing the return code of cwltool, not the tool it is running, so we don't need to
    // worry about success/temporaryFail/permanentFail codes.
    (true, Some(Set(0)))
  }

  private lazy val outputParams: Map[DxName, OutputParameter] = {
    val outputParams: Map[DxName, OutputParameter] = tool.outputs.flatMap { param =>
      param.id.map(_.name).map(CwlDxName.fromSourceName(_) -> param)
    }.toMap
    if (logger.isVerbose) {
      logger.traceLimited(s"outputParams=${outputParams}")
    }
    outputParams
  }

  override protected def evaluateOutputs(
      localizedInputs: Map[DxName, (Type, Value)]
  ): (Map[DxName, (Type, Value)], Map[DxName, (Set[String], Map[String, String])]) = {
    // the outputs were written to stdout
    val stdoutFile = workerPaths.getStdoutFile().asJavaPath
    val localizedOutputs = if (Files.exists(stdoutFile)) {
      val allOutputs: Map[DxName, JsValue] = JsUtils.jsFromFile(stdoutFile) match {
        case JsObject(outputs) =>
          outputs.map {
            case (name, jsv) => CwlDxName.fromDecodedName(name) -> jsv
          }
        case JsNull => Map.empty[DxName, JsValue]
        case other  => throw new Exception(s"unexpected cwltool outputs ${other}")
      }
      if (logger.isVerbose) {
        trace(s"cwltool outputs:\n${JsObject(allOutputs.map {
          case (dxName, jsv) => dxName.decoded -> jsv
        }).prettyPrint}")
      }
      val cwlOutputs = outputParams.map {
        case (name, param) if allOutputs.contains(name) =>
          name -> CwlValue.deserialize(allOutputs(name), param.cwlType, typeAliases)
        case (name, param) if CwlOptional.isOptional(param.cwlType) =>
          name -> (param.cwlType, NullValue)
        case (_, param) =>
          throw new Exception(s"missing value for output parameter ${param}")
      }
      CwlUtils.toIR(cwlOutputs)
    } else {
      Map.empty[DxName, (Type, Value)]
    }
    (localizedOutputs, Map.empty)
  }

  override protected def outputTypes: Map[DxName, Type] = {
    outputParams.map {
      case (name, param) => name -> CwlUtils.toIRType(param.cwlType)
    }
  }
}
