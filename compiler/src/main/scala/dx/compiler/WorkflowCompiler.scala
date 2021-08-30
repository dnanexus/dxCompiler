package dx.compiler

import dx.api.{DxApi, DxUtils, DxWorkflowStage, Field}
import dx.core.Constants
import dx.core.ir._
import dx.translator.CallableAttributes._
import dx.translator.Extras
import dx.util.{CodecUtils, FileSourceResolver, Logger}
import spray.json._

import scala.collection.immutable.TreeSeqMap
import dx.translator.DockerRegistry
import dx.api.DxWorkflow
import dx.api.DxApplet

case class WorkflowCompiler(separateOutputs: Boolean,
                            extras: Option[Extras],
                            parameterLinkSerializer: ParameterLinkSerializer,
                            useManifests: Boolean,
                            complexPathValues: Boolean,
                            fileResolver: FileSourceResolver,
                            dxApi: DxApi = DxApi.get,
                            logger: Logger = Logger.get)
    extends ExecutableCompiler(extras,
                               parameterLinkSerializer,
                               complexPathValues,
                               fileResolver,
                               dxApi) {

  private def workflowInputParameterToNative(parameter: Parameter,
                                             stageInput: StageInput): Vector[JsValue] = {
    // The default value can come from the stage input or the workflow input
    def getDefault(stageInput: StageInput): Option[Option[Value]] = {
      stageInput match {
        case WorkflowInput(wfParam) if wfParam != parameter => Some(wfParam.defaultValue)
        case StaticInput(value)                             => Some(Some(value))
        case ArrayInput(inputs) =>
          val inputValues = inputs.flatMap(getDefault)
          if (inputValues.isEmpty) {
            None
          } else {
            Some(Some(Value.VArray(inputValues.flatten)))
          }
        case _ => None
      }
    }
    val paramWithDefault = getDefault(stageInput)
      .map(value => parameter.copy(defaultValue = value))
      .getOrElse(parameter)
    inputParameterToNative(paramWithDefault)
  }

  /**
    * Converts a workflow output paramter to outputSpec JSON. A complex output
    * will generate two output specs - one for the value and one for an array
    * of files that are nested within the complex value.
    * @param parameter the workflow output Parameter
    * @param stageInput the StageInput from which the workflow output is derived
    * @return
    */
  private def workflowOutputParameterToNative(parameter: Parameter,
                                              stageInput: StageInput): Vector[JsValue] = {
    val outputSpec: Map[String, Map[String, JsValue]] =
      outputParameterToNative(parameter).map { obj =>
        val name = obj.fields.get("name") match {
          case Some(JsString(name)) => name
          case other =>
            throw new Exception(s"Unexpected value for 'name' field: ${other}")
        }
        name -> obj.fields
      }.toMap

    def createLink(stageInput: StageInput, dxType: Type): ParameterLink = {
      stageInput match {
        case StaticInput(value) => parameterLinkSerializer.createLink(dxType, value)
        case LinkInput(dxStage, paramName) =>
          ParameterLinkStage(dxStage, IORef.Output, paramName, dxType)
        case WorkflowInput(wfParam) =>
          // TODO: if the input has a non-static default, link to the value of the workflow input
          //  (either the user-specified value or the result of evaluting the expression)
          //  - right now this only links to the user-specified value
          ParameterLinkWorkflowInput(wfParam.name, dxType)
        case ArrayInput(stageInputs) =>
          val itemType = dxType match {
            case Type.TArray(itemType, _) if Type.isNative(itemType, !complexPathValues) => itemType
            case _ =>
              throw new Exception(
                  s"""ArrayInput for workflow output parameter ${parameter.name} with type ${dxType} 
                     |that is not a natively supported array type""".stripMargin
                    .replaceAll("\n", " ")
              )
          }
          ParameterLinkValue(JsArray(stageInputs.map { i =>
            parameterLinkSerializer.serializeSimpleLink(createLink(i, itemType))
          }), dxType)
        case other =>
          throw new Exception(s"Bad value for stage input ${other}")
      }
    }

    val outputSources: Vector[(DxName, JsValue)] =
      parameterLinkSerializer.createFieldsFromLink(createLink(stageInput, parameter.dxType),
                                                   parameter.name)

    // merge the specification and the output sources
    outputSources.map {
      case (fieldName, output) =>
        JsObject(outputSpec(fieldName.encoded) ++ Map("outputSource" -> output))
    }
  }

  /**
    * Translates stage inputs to native input spec.
    * @param inputs stage inputs
    * @return
    */
  private def stageInputToNative(inputs: Vector[(Parameter, StageInput)]): JsValue = {
    def createLink(stageInput: StageInput, dxType: Type): Option[ParameterLink] = {
      stageInput match {
        case EmptyInput =>
          // We do not have a value for this input at compile time. For compulsory applet inputs,
          // the user will have to fill in a value at runtime.
          None
        case StaticInput(value) =>
          Some(parameterLinkSerializer.createLink(dxType, value))
        case LinkInput(dxStage, paramName) =>
          Some(ParameterLinkStage(dxStage, IORef.Output, paramName, dxType))
        case WorkflowInput(wfParam) =>
          Some(ParameterLinkWorkflowInput(wfParam.name, dxType))
        case ArrayInput(stageInputs) =>
          val itemType = dxType match {
            case Type.TArray(itemType, _) if Type.isNative(itemType, !complexPathValues) => itemType
            case _ =>
              throw new Exception(
                  s"ArrayInput for stage input with type ${dxType} that is not a natively supported array type"
              )
          }
          val inputs = stageInputs.flatMap { i =>
            createLink(i, itemType).map(parameterLinkSerializer.serializeSimpleLink)
          }
          if (inputs.isEmpty) {
            None
          } else {
            Some(ParameterLinkValue(JsArray(inputs), dxType))
          }
      }
    }
    // sort the inputs, to make the request deterministic
    JsObject(
        inputs
          .flatMap {
            case (parameter, stageInput) =>
              createLink(stageInput, parameter.dxType)
                .map { link =>
                  parameterLinkSerializer.createFieldsFromLink(link, parameter.name)
                }
                .getOrElse(Vector.empty)
          }
          .map {
            case (dxName, jsv) => dxName.encoded -> jsv
          }
          .to(TreeSeqMap)
    )
  }

  private def workflowAttributesToNative(
      workflow: Workflow,
      defaultTags: Set[String],
      extendedDescription: Option[String]
  ): (Map[String, JsValue], Map[String, JsValue]) = {
    val (commonMeta, commonDetails) =
      callableAttributesToNative(
          callable = workflow,
          defaultTags = defaultTags,
          extendedDescription = extendedDescription
      )
    val workflowMeta = workflow.attributes.collect {
      // These will be implemented in a future PR
      case _: CallNamesAttribute =>
        throw new NotImplementedError()
      case _: RunOnSingleNodeAttribute =>
        throw new NotImplementedError()
      // These are currently ignored because they only apply to apps
      //case VersionAttribute(text) => Some("version" -> JsString(text))
    }
    (commonMeta ++ workflowMeta, commonDetails)
  }

  // Summarize reportable dependencies from details of a workflow
  private def reportableDependenciesFromWorkflow(
      workflow: DxWorkflow
  ): Option[String] = {
    val description = workflow.describe(Set(Field.Name, Field.Details, Field.Hidden))

    // Do not generate description section if hidden
    if (description.hidden.exists(h => h)) {
      return None
    }

    // Generate Markdown description section
    val headerMd = s"\n## Sub-workflow ${description.name}"
    val md = description.details match {
      case Some(JsObject(d)) => {
        val appsMd = d.get(Constants.NativeAppDependencies) match {
          case Some(JsArray(a)) =>
            s"${listMd("Native apps")}${a.map(v => listMd2(v.toString)).mkString}"
          case _ => ""
        }
        val filesMd = d.get(Constants.FileDependencies) match {
          case Some(JsArray(a)) =>
            s"${listMd("Files")}${a.map(v => listMd2(v.toString)).mkString}"
          case _ => ""
        }
        val dockerRegistryCredentialsMd =
          d.get(Constants.DockerRegistryCredentialsUri) match {
            case Some(JsString(s)) => s"${listMd("Docker registry credentials file")}${listMd2(s)}"
            case _                 => ""
          }
        s"${appsMd}${filesMd}${dockerRegistryCredentialsMd}"
      }
      case _ => ""
    }
    Option.when(md.nonEmpty)(s"${headerMd}${md}")
  }

  // Summarize reportable dependencies from details of an applet
  private def reportableDependenciesFromApplication(
      applet: DxApplet
  ): Option[String] = {
    val description = applet.describe(Set(Field.Name, Field.Details, Field.Hidden))

    // Do not generate description section if hidden
    if (description.hidden.exists(h => h)) {
      return None
    }

    // Generate Markdown description section
    val headerMd = s"\n## Task ${description.name}"
    val md = description.details match {
      case Some(JsObject(d)) => {
        val filesMd = d.get(Constants.FileDependencies) match {
          case Some(JsArray(a)) =>
            s"${listMd("Files")}${a.map(v => listMd2(v.toString)).mkString}"
          case _ => ""
        }
        val networkDockerImageMd = d.get(Constants.NetworkDockerImage) match {
          case Some(JsString(s)) => s"${listMd("Network Docker image")}${listMd2(s)}"
          case _                 => ""
        }
        val dynamicDockerImageMd = d.get(Constants.DynamicDockerImage) match {
          case Some(JsBoolean(b)) if b => listMd("Docker image determined at runtime")
          case _                       => ""
        }
        val staticInstanceTypeMd = d.get(Constants.StaticInstanceType) match {
          case Some(JsString(s)) => s"${listMd("Hard-coded instance type")}${listMd2(s)}"
          case _                 => ""
        }
        s"${filesMd}${networkDockerImageMd}${dynamicDockerImageMd}${staticInstanceTypeMd}"
      }
      case _ => ""
    }
    Option.when(md.nonEmpty)(s"${headerMd}${md}")
  }

  // TODO: Use templates for Markdown dependency report

  // Summarize dependencies that are not bundled for cloning
  // with the workflow for adding to description in Markdown
  private def summarizeReportableDependencies(
      workflow: Workflow,
      executableDict: Map[String, CompiledExecutable],
      execTree: JsObject,
      nativeAppDependencies: Vector[String],
      fileDependencies: Vector[String],
      dockerRegistryCredentialsUri: Option[String]
  ): Option[String] = {
    val headerMd =
      "# This workflow requires access to the following dependencies that are not " +
        "packaged with the workflow"

    val topLevelMd = Vector(
        Option.when(nativeAppDependencies.nonEmpty)(
            Vector(listMd("Native apps")) ++ nativeAppDependencies.map(listMd2)
        ),
        Option.when(fileDependencies.nonEmpty)(
            Vector(listMd("Files")) ++ fileDependencies.map(listMd2)
        ),
        Option.when(dockerRegistryCredentialsUri.nonEmpty)(
            Vector(listMd("Docker registry credentials file"),
                   listMd2(dockerRegistryCredentialsUri.get))
        )
    ).flatten.flatten.mkString

    val executableNames = ExecutableTree.extractExecutableNames(execTree)
    val executables = executableNames.collect(executableDict.get(_)).flatten

    val subWorkflowsMd = executables
      .collect {
        case CompiledExecutable(_: Workflow, dxExec: DxWorkflow, _, _) =>
          dxExec
      }
      .flatMap(reportableDependenciesFromWorkflow)
      .mkString

    val appletsMd = executables
      .collect {
        case CompiledExecutable(Application(_, _, _, _, _, ExecutableKindApplet, _, _, _, _, _, _),
                                dxExec: DxApplet,
                                _,
                                _) =>
          dxExec
      }
      .flatMap(reportableDependenciesFromApplication)
      .mkString

    // Append header & return full report if it has contents
    val md = s"${topLevelMd}${subWorkflowsMd}${appletsMd}"
    Option.when(md.nonEmpty)(s"${headerMd}${md}")
  }

  def apply(
      workflow: Workflow,
      executableDict: Map[String, CompiledExecutable]
  ): (Map[String, JsValue], JsValue) = {
    logger.trace(s"Building /workflow/new request for ${workflow.name}")

    // build the "inputs" and "outputs" part of the API request
    val wfInputOutput: Map[String, JsValue] = {
      // Locked workflows have well defined inputs and outputs
      if (workflow.locked) {
        val inputParams = if (useManifests) {
          // When using manifests, there are two input parameters:
          // * "input_manifests___": an array of all the manifest files from the
          //   upstream stages in any workflow LinkInputs
          // * "input_links___": a hash of parameter names to values, where a value
          //   may be a static value, which must be embedded in a hash with a "value___"
          //   key, or a link to a value in one of the manifests, which is a hash with
          //   a single field mapping the stage ID to the value. A value can also be an
          //   array of nested value/link hashes.
          def getWorkflowInputValue(stageInput: StageInput,
                                    dxType: Type): (Set[DxWorkflowStage], Option[Value]) = {
            stageInput match {
              case EmptyInput => (Set.empty, None)
              case StaticInput(value) =>
                (Set.empty, Some(Value.VHash(Constants.ValueKey.decoded -> value)))
              case LinkInput(stageId, stageParamName) =>
                (Set(stageId),
                 Some(Value.VHash(stageId.id -> Value.VString(stageParamName.decoded))))
              case ArrayInput(inputs) =>
                val itemType = dxType match {
                  case Type.TArray(itemType, _) => itemType
                  case _ =>
                    throw new Exception(
                        s"ArrayInput for stage input with type ${dxType} that is not a natively supported array type"
                    )
                }
                val (stages, values) = inputs.map(getWorkflowInputValue(_, itemType)).unzip
                (stages.flatten.toSet, Some(Value.VArray(values.flatten)))
              case _ =>
                throw new Exception(s"invalid nested stage input ${stageInput}")
            }
          }
          val (inputStages, inputLinks) = workflow.inputs.map {
            case (param, EmptyInput)
                if workflow.level == Level.Top ||
                  param.defaultValue.isDefined ||
                  Type.isOptional(param.dxType) =>
              // this is one of:
              // 1) an optional input - the value will come from the default or will be null
              // 2) a top-level workflow input - the value (if any) will be supplied in the
              //    input manifest
              (Set.empty, None)
            case (param, EmptyInput) =>
              // empty required inputs are not allowed for nested locked workflows
              throw new Exception(
                  s"locked workflow ${workflow.name} has empty required input ${param.name}"
              )
            case (param, WorkflowInput(wfParam)) if param == wfParam =>
              // this is one of:
              // 1) an optional input - the value will come from the default or will be null
              // 2) a top-level workflow input - the value (if any) will be supplied in the
              //    input manifest
              // 3) a scatter variable (or an expression that references the scatter
              //    variable) - the value is provided by the caller
              (Set.empty, None)
            case (param, stageInput) =>
              val (inputStages, value) = getWorkflowInputValue(stageInput, param.dxType)
              (inputStages, value.map(v => param.name.decoded -> v))
          }.unzip
          val stageManifestLinks: StageInput = ArrayInput(inputStages.flatten.toSet.map { stage =>
            LinkInput(stage, Constants.OutputManifest)
          }.toVector)
          Vector(
              (ExecutableCompiler.InputManifestParameter, EmptyInput),
              (ExecutableCompiler.InputManfestFilesParameter, stageManifestLinks),
              (ExecutableCompiler.InputLinksParameter,
               StaticInput(Value.VHash(inputLinks.flatten.to(TreeSeqMap)))),
              (ExecutableCompiler.OutputIdParameter, StaticInput(Value.VString(workflow.name))),
              (ExecutableCompiler.CallNameParameter, EmptyInput)
          )
        } else {
          workflow.inputs
        }

        val outputParams = if (useManifests) {
          // When using manifests, there is a single "output_manifest___" parameter, which
          // is a manifest file. The common output applet is responsible for merging all
          // incoming manifests into a single manifest, so the workflow output is just a
          // link to the common stage output.
          val outputStage = workflow.stages.filter(_.description == Constants.OutputStage) match {
            case Vector(outputStage) => outputStage
            case _                   => throw new Exception("expected a single output stage")
          }
          Vector(
              (ExecutableCompiler.OutputManifestParameter,
               LinkInput(outputStage.dxStage, Constants.OutputManifest))
          )
        } else {
          workflow.outputs
        }

        val wfInputSpec: Vector[JsValue] = inputParams
          .sortWith(_._1.name < _._1.name)
          .flatMap {
            case (parameter, stageInput) => workflowInputParameterToNative(parameter, stageInput)
          }
        val wfOutputSpec: Vector[JsValue] = outputParams
          .sortWith(_._1.name < _._1.name)
          .flatMap {
            case (parameter, stageInput) => workflowOutputParameterToNative(parameter, stageInput)
          }
        Map("inputs" -> JsArray(wfInputSpec), "outputs" -> JsArray(wfOutputSpec))
      } else {
        Map.empty
      }
    }

    // build the "stages" part of the API request
    val stages =
      workflow.stages.map { stage =>
        val CompiledExecutable(irExecutable, dxExec, _, _) = executableDict(stage.calleeName)
        val linkedInputs = if (useManifests) {
          // when using manifests, we have to create an input array of all the
          // manifests output by any linked stages, and a hash of links between
          // the manifest fields and the stage inputs, plus any static values
          def getStageInputValue(
              stageInput: StageInput,
              dxType: Type
          ): (Boolean, Set[DxWorkflowStage], Option[Value]) = {
            stageInput match {
              case EmptyInput => (false, Set.empty, None)
              case StaticInput(value) =>
                (false, Set.empty, Some(Value.VHash(Constants.ValueKey.decoded -> value)))
              case WorkflowInput(wfParam) =>
                (true,
                 Set.empty,
                 Some(
                     Value.VHash(
                         Constants.WorkflowKey.decoded -> Value.VString(wfParam.name.decoded)
                     )
                 ))
              case LinkInput(sourceStage, sourceParamName) =>
                (false,
                 Set(sourceStage),
                 Some(Value.VHash(sourceStage.id -> Value.VString(sourceParamName.decoded))))
              case ArrayInput(inputs) =>
                val itemType = dxType match {
                  case Type.TArray(itemType, _) => itemType
                  case _ =>
                    throw new Exception(
                        s"ArrayInput for stage input with type ${dxType} that is not a natively supported array type"
                    )
                }
                val (inputWorkflow, inputStages, values) =
                  inputs.map(getStageInputValue(_, itemType)).unzip3
                (inputWorkflow.exists(identity),
                 inputStages.flatten.toSet,
                 Some(Value.VArray(values.flatten)))
            }
          }
          val (inputWorkflow, inputStages, inputLinks) = irExecutable.inputVars
            .zip(stage.inputs)
            .map {
              case (param, EmptyInput)
                  if Type.isOptional(param.dxType) || param.defaultValue.isDefined =>
                (false, Set.empty, None)
              case (param, EmptyInput) =>
                throw new Exception(
                    s"""no input specified for stage input ${stage}.${param.name}; overriding required
                       |call inputs at runtime is not supported when using manifests""".stripMargin
                      .replaceAll("\n", " ")
                )
              case (param, stageInput) =>
                val (inputWorkflow, inputStages, value) =
                  getStageInputValue(stageInput, param.dxType)
                (inputWorkflow, inputStages, value.map(v => param.name.decoded -> v))
            }
            .unzip3
          val inputStageManifests = inputStages.flatten.toSet.map { stage =>
            LinkInput(stage, Constants.OutputManifest)
          }.toVector
          // the manifest ID for the output stage comes from the workflow
          val outputInputs = if (stage.description == Constants.OutputStage) {
            Vector(
                (ExecutableCompiler.OutputIdParameter,
                 WorkflowInput(ExecutableCompiler.OutputIdParameter)),
                (ExecutableCompiler.CallNameParameter,
                 WorkflowInput(ExecutableCompiler.CallNameParameter))
            )
          } else {
            Vector(
                (ExecutableCompiler.OutputIdParameter, StaticInput(Value.VString(stage.dxStage.id)))
            )
          }
          val stageInputs = Vector(
              (ExecutableCompiler.InputManfestFilesParameter, ArrayInput(inputStageManifests)),
              (ExecutableCompiler.InputLinksParameter,
               StaticInput(Value.VHash(inputLinks.flatten.to(TreeSeqMap))))
          ) ++ outputInputs
          // if there are any workflow inputs, pass the workflow manifests and links to the stage
          val stageWorkflowInputs = if (inputWorkflow.exists(identity)) {
            Vector(
                (ExecutableCompiler.WorkflowInputManifestParameter,
                 WorkflowInput(ExecutableCompiler.InputManifestParameter)),
                (ExecutableCompiler.WorkflowInputManfestFilesParameter,
                 WorkflowInput(ExecutableCompiler.InputManfestFilesParameter)),
                (ExecutableCompiler.WorkflowInputLinksParameter,
                 WorkflowInput(ExecutableCompiler.InputLinksParameter))
            )
          } else {
            Vector.empty
          }
          stageInputs ++ stageWorkflowInputs
        } else {
          irExecutable.inputVars.zip(stage.inputs)
        }
        val stageAttrs = Vector(
            "id" -> JsString(stage.dxStage.id),
            "executable" -> JsString(dxExec.id),
            "name" -> JsString(stage.description),
            "input" -> stageInputToNative(linkedInputs)
        )
        val folderAttr = if (separateOutputs) {
          val folderName = s"${workflow.name}_${stage.dxStage.id}"
          Vector("folder" -> JsString(folderName))
        } else {
          Vector.empty
        }
        // convert the per-stage metadata into JSON
        JsObject((stageAttrs ++ folderAttr).toMap)
      }
    // build the details JSON
    val defaultTags = Set(Constants.CompilerTag)
    // compress and base64 encode the source code
    val sourceDetails = Map(
        Constants.SourceCode -> JsString(
            CodecUtils.gzipAndBase64Encode(workflow.document.toString)
        ),
        Constants.ParseOptions -> workflow.document.optionsToJson
    )
    // links through applets that run workflow fragments
    val transitiveDependencies: Vector[ExecutableLink] =
      workflow.stages.flatMap { stage =>
        val CompiledExecutable(_, _, dependencies, _) = executableDict(stage.calleeName)
        dependencies
      }
    // link to applets used by the fragments. This notifies the platform that they
    // need to be cloned when copying workflows.
    val dxLinks = transitiveDependencies.map {
      case ExecutableLink(name, _, _, dxExec) =>
        s"link_${name}" -> JsObject(DxUtils.DxLinkKey -> JsString(dxExec.id))
    }.toMap
    val delayDetails = delayWorkspaceDestructionToNative
    // generate the executable tree
    val execTree = ExecutableTree(executableDict).apply(workflow)
    val execTreeDetails = Map(
        "execTree" -> JsString(CodecUtils.gzipAndBase64Encode(execTree.toString))
    )

    // Collect native app dependencies
    val nativeAppDependencies: Vector[String] =
      Set
        .from(
            workflow.stages
              .collect { stage =>
                executableDict(stage.calleeName) match {
                  case CompiledExecutable(
                      Application(_,
                                  _,
                                  _,
                                  _,
                                  _,
                                  ExecutableKindNative(_, _, _, _, _),
                                  _,
                                  _,
                                  _,
                                  _,
                                  _,
                                  _),
                      dxExec,
                      _,
                      _
                      ) =>
                    dxExec.id
                }
              }
              .filter(_.startsWith("app-"))
        )
        .toVector
    val nativeAppDependenciesDetails = Option
      .when(nativeAppDependencies.nonEmpty)(
          Map(
              Constants.NativeAppDependencies -> JsArray(
                  nativeAppDependencies.map(id => JsString(id))
              )
          )
      )
      .getOrElse(Map.empty)

    // Collect file dependencies from inputs & private variables
    val fileInputParams = workflow.inputs
      .map(i => i._1)
      .filter(param => Type.unwrapOptional(param.dxType) == Type.TFile)
    val fileDependencies = (fileInputParams
      .flatMap(fileDependenciesFromParam) ++ workflow.staticFileDependencies).distinct
    val fileDependenciesDetails = Option
      .when(fileDependencies.nonEmpty)(
          Map(Constants.FileDependencies -> JsArray(fileDependencies.map(s => JsString(s))))
      )
      .getOrElse(Map.empty)

    // Collect Docker registry credentials URI
    val dockerRegistryCredentialsUri: Option[String] = extras.flatMap(_.dockerRegistry) match {
      case Some(DockerRegistry(_, credentials, _, _)) => Some(credentials)
      case _                                          => None
    }
    val dockerRegistryCredentialsUriDetails = dockerRegistryCredentialsUri
      .map(cred => Map(Constants.DockerRegistryCredentialsUri -> JsString(cred)))
      .getOrElse(Map.empty)

    // Build extended description with dependency report
    val extendedDescription =
      summarizeReportableDependencies(
          workflow = workflow,
          executableDict = executableDict,
          execTree = execTree.asJsObject,
          nativeAppDependencies = nativeAppDependencies,
          fileDependencies = fileDependencies,
          dockerRegistryCredentialsUri = dockerRegistryCredentialsUri
      )
    val (wfMeta, wfMetaDetails) =
      workflowAttributesToNative(workflow, defaultTags, extendedDescription)

    // Collect details
    val details = wfMetaDetails ++ sourceDetails ++ dxLinks ++ delayDetails ++ execTreeDetails ++
      nativeAppDependenciesDetails ++ fileDependenciesDetails ++ dockerRegistryCredentialsUriDetails
    val isTopLevel = workflow.level == Level.Top
    // required arguments for API call
    val required = Map(
        "name" -> JsString(workflow.name),
        "stages" -> JsArray(stages),
        "details" -> JsObject(details),
        // Sub-workflow are compiled to hidden objects.
        "hidden" -> JsBoolean(!isTopLevel)
    )
    // look for ignoreReuse in runtime hints and in extras - the later overrides the former
    val ignoreReuse = extras.flatMap(_.ignoreReuse) match {
      case Some(true) => Map("ignoreReuse" -> JsArray(JsString("*")))
      case _          => Map.empty
    }
    // build API request
    val request = required ++ wfInputOutput ++ wfMeta ++ ignoreReuse
    (request, execTree)
  }
}
