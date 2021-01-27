package dx.compiler

import dx.api.{DxApi, DxUtils, DxWorkflowStage}
import dx.core.Constants
import dx.core.ir._
import dx.translator.CallableAttributes._
import dx.translator.Extras
import dx.util.{CodecUtils, Logger}
import spray.json._

import scala.collection.immutable.{SeqMap, TreeMap, TreeSeqMap}

case class WorkflowCompiler(extras: Option[Extras],
                            parameterLinkSerializer: ParameterLinkSerializer,
                            useManifests: Boolean,
                            dxApi: DxApi = DxApi.get,
                            logger: Logger = Logger.get)
    extends ExecutableCompiler(extras, parameterLinkSerializer, dxApi) {

  private def workflowInputParameterToNative(parameter: Parameter,
                                             stageInput: StageInput): Vector[JsValue] = {
    // The default value can come from the stageInput or the parameter
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
    * Converts a workflow output paramter to outputSpec JSON. A complex ouptut
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
          ParameterLinkWorkflowInput(wfParam.dxName, dxType)
        case ArrayInput(stageInputs) =>
          val itemType = dxType match {
            case Type.TArray(itemType, _) if Type.isNative(itemType) => itemType
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

    val outputSources: Vector[(String, JsValue)] =
      parameterLinkSerializer.createFieldsFromLink(createLink(stageInput, parameter.dxType),
                                                   parameter.dxName)

    // merge the specification and the output sources
    outputSources.map {
      case (fieldName, output) =>
        JsObject(outputSpec(fieldName) ++ Map("outputSource" -> output))
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
          Some(ParameterLinkWorkflowInput(wfParam.dxName, dxType))
        case ArrayInput(stageInputs) =>
          val itemType = dxType match {
            case Type.TArray(itemType, _) if Type.isNative(itemType) => itemType
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
    val jsInputs: TreeMap[String, JsValue] = inputs.foldLeft(TreeMap.empty[String, JsValue]) {
      case (accu, (parameter, stageInput)) =>
        accu ++ createLink(stageInput, parameter.dxType)
          .map { link =>
            parameterLinkSerializer.createFieldsFromLink(link, parameter.dxName).to(TreeSeqMap)
          }
          .getOrElse(SeqMap.empty)
    }
    JsObject(jsInputs)
  }

  private def workflowAttributesToNative(
      workflow: Workflow,
      defaultTags: Set[String]
  ): (Map[String, JsValue], Map[String, JsValue]) = {
    val (commonMeta, commonDetails) = callableAttributesToNative(workflow, defaultTags)
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

  def apply(
      workflow: Workflow,
      executableDict: Map[String, CompiledExecutable]
  ): (Map[String, JsValue], JsValue) = {
    logger.trace(s"Building /workflow/new request for ${workflow.name}")

    // convert inputs and outputs to spec
    val wfInputOutput: Map[String, JsValue] = {
      if (workflow.locked) {
        // Locked workflows have well defined inputs and outputs
        val (inputParams, outputParams) = if (useManifests) {
          // When using manifests, there are two input parameters:
          // * "input_manifests___": an array of all the manifest files from the
          //   upstream stages in any workflow LinkInputs
          // * "input_links___": a hash of parameter names to values, where a value
          //   may be a static value, which must be embedded in a hash with a "value___"
          //   key, or a link to a value in one of the manifests, which is a hash with
          //   a single field mapping the stage ID to the value. A value can also be an
          //   array of nested value/link hashes.
          def getStageInputValue(stageInput: StageInput,
                                 dxType: Type): (Set[DxWorkflowStage], Option[Value]) = {
            stageInput match {
              case EmptyInput => (Set.empty, None)
              case StaticInput(value) =>
                (Set.empty, Some(Value.VHash(TreeSeqMap(Constants.ValueKey -> value))))
              case LinkInput(stageId, stageParamName) =>
                (Set(stageId),
                 Some(Value.VHash(TreeSeqMap(stageId.id -> Value.VString(stageParamName)))))
              case ArrayInput(inputs) =>
                val itemType = dxType match {
                  case Type.TArray(itemType, _) => itemType
                  case _ =>
                    throw new Exception(
                        s"ArrayInput for stage input with type ${dxType} that is not a natively supported array type"
                    )
                }
                val (stages, values) = inputs.map(getStageInputValue(_, itemType)).unzip
                (stages.flatten.toSet, Some(Value.VArray(values.flatten)))
              case _ =>
                throw new Exception(s"invalid nested stage input ${stageInput}")
            }
          }
          val (inputStages, inputLinks) = workflow.inputs.map {
            case (param, EmptyInput)
                if param.defaultValue.isDefined || Type.isOptional(param.dxType) =>
              // an optional input - the value will come from the default or will be null
              (Set.empty, None)
            case (_, EmptyInput) if workflow.level == Level.Top =>
              // a top-level workflow input - the value (if any) will be supplied
              // in the input manifest
              (Set.empty, None)
            case (param, WorkflowInput(wfParam))
                if param == wfParam && (
                    param.defaultValue.isDefined || Type.isOptional(param.dxType)
                ) =>
              // an optional input - the value will come from the default or will be null
              (Set.empty, None)
            case (param, WorkflowInput(wfParam))
                if param == wfParam && workflow.level == Level.Top =>
              // a top-level workflow input - this is a top-level workflow,
              // so the value (if any) will be supplied in the input manifest
              (Set.empty, None)
            case (param, EmptyInput | _: WorkflowInput) =>
              // empty required inputs are not allowed for nested locked workflows
              throw new Exception(
                  s"locked workflow ${workflow.name} has empty required input ${param.name}"
              )
            case (param, stageInput) =>
              val (inputStages, value) = getStageInputValue(stageInput, param.dxType)
              (inputStages, value.map(v => param.name -> v))
          }.unzip
          val inputParams = Vector(
              (Parameter(Constants.InputManifests, Type.TArray(Type.TFile)),
               ArrayInput(inputStages.flatten.toSet.map { stage =>
                 LinkInput(stage, Constants.OutputManifest)
               }.toVector)),
              (Parameter(Constants.InputLinks, Type.THash),
               StaticInput(Value.VHash(inputLinks.flatten.to(TreeSeqMap))))
          )
          // When using manifests, there is a single "output_manifest___" parameter, which
          // is a manifest file. The common output applet is responsible for merging all
          // incomming manifests into a single manifest, so the workflow output is just a
          // link to the common stage output.
          val outputStage = workflow.stages.filter(_.description == Constants.OutputStage) match {
            case Vector(outputStage) => outputStage
            case _                   => throw new Exception("expected a single output stage")
          }
          val outputParams = Vector(
              (Parameter(Constants.OutputManifest, Type.TFile),
               LinkInput(outputStage.dxStage, Constants.OutputManifest))
          )
          (inputParams, outputParams)
        } else {
          (workflow.inputs, workflow.outputs)
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

    val stages =
      workflow.stages.foldLeft(Vector.empty[JsValue]) {
        case (stagesReq, stg) =>
          val CompiledExecutable(irApplet, dxExec, _, _) = executableDict(stg.calleeName)
          val linkedInputs = if (useManifests) {
            // when using manifests, we have to create an input array of all the
            // manifests output by any linked stages, and a hash of links between
            // the manifest fields and the stage inputs, plus any static values
            val (wfLinks, stageLinks, manifests, staticInputs) = irApplet.inputVars
              .zip(stg.inputs)
              .foldLeft(Map.empty[String, JsValue],
                        Map.empty[String, JsValue],
                        Set.empty[x],
                        Vector.empty[x]) {
                case (accu, (param, EmptyInput))
                    if Type.isOptional(param.dxType) || param.defaultValue.isDefined =>
                  accu
                case (_, (param, EmptyInput)) =>
                  throw new Exception(
                      s"""no input specified for stage input ${stg}.${param.name}; overriding required 
                         |call inputs at runtime is not supported when using manifests""".stripMargin
                        .replaceAll("\n", " ")
                  )
                case ((wfAccu, stageAccu, manifestAccu, staticAccu),
                      (param, wfInput: WorkflowInput)) =>
                  (wfAccu + (param.name -> JsString(wfInput.param.name)),
                   stageAccu,
                   manifestAccu +,
                   staticAccu)
                case ((wfAccu, linkAccu, staticAccu), (param, wfInput: LinkInput)) =>
                  (wfAccu, linkAccu :+ (param, wfInput), staticAccu)
                case ((wfAccu, linkAccu, staticAccu), (param, staticInput: StaticInput)) =>
                  (wfAccu, linkAccu, staticAccu :+ (param, staticInput))
              }
          } else {
            irApplet.inputVars.zip(stg.inputs)
          }
          // convert the per-stage metadata into JSON
          val stageReqDesc = JsObject(
              Map("id" -> JsString(stg.dxStage.id),
                  "executable" -> JsString(dxExec.id),
                  "name" -> JsString(stg.description),
                  "input" -> stageInputToNative(linkedInputs))
          )
          stagesReq :+ stageReqDesc
      }
    // build the details JSON
    val defaultTags = Set(Constants.CompilerTag)
    val (wfMeta, wfMetaDetails) = workflowAttributesToNative(workflow, defaultTags)
    // compress and base64 encode the source code
    val sourceDetails = Map(
        Constants.SourceCode -> JsString(
            CodecUtils.gzipAndBase64Encode(workflow.document.toString)
        )
    )
    // links through applets that run workflow fragments
    val transitiveDependencies: Vector[ExecutableLink] =
      workflow.stages.foldLeft(Vector.empty[ExecutableLink]) {
        case (accu, stg) =>
          val CompiledExecutable(_, _, dependencies, _) = executableDict(stg.calleeName)
          accu ++ dependencies
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
    val details = wfMetaDetails ++ sourceDetails ++ dxLinks ++ delayDetails ++ execTreeDetails
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
