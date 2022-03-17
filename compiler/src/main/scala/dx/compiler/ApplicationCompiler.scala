package dx.compiler

import dx.api.{
  DxAccessLevel,
  DxApi,
  DxApplet,
  DxFileDescribe,
  DxPath,
  DxProject,
  DxUtils,
  DxWorkflow,
  InstanceTypeRequest
}
import dx.core.Constants
import dx.core.io.{DxWorkerPaths, StreamFiles}
import dx.core.ir._
import dx.core.ir.RunSpec._
import dx.translator.{DockerRegistry, DxAccess, DxRunSpec, DxTimeout, Extras}
import dx.translator.CallableAttributes._
import dx.translator.ExtrasJsonProtocol._
import dx.util.{CodecUtils, FileSourceResolver, Logger}
import spray.json._
import wdlTools.generators.Renderer

import scala.util.{Failure, Success, Try}

object ApplicationCompiler {
  val DefaultAppletTimeoutInDays = 2
  // templates
  private val GenericDockerPreambleTemplate = "templates/generic_docker_preamble.ssp"
  private val EcrDockerPreambleTemplate = "templates/ecr_docker_preamble.ssp"
  private val DynamicAppletCommandTemplate = "templates/dynamic_applet_script.ssp"
  private val StaticAppletCommandTemplate = "templates/static_applet_script.ssp"
  private val WorkflowFragmentTemplate = "templates/workflow_fragment_script.ssp"
  private val WorkflowCommandTemplate = "templates/workflow_command_script.ssp"
  // keys used in templates
  private val RegistryKey = "registry"
  private val CredentialsKey = "credentials"
  private val UsernameKey = "username"
  private val AwsRegionKey = "region"
}

case class ApplicationCompiler(typeAliases: Map[String, Type],
                               runtimeAsset: Option[JsValue],
                               runtimeJar: String,
                               runtimePathConfig: DxWorkerPaths,
                               runtimeTraceLevel: Int,
                               separateOutputs: Boolean,
                               streamFiles: StreamFiles.StreamFiles,
                               waitOnUpload: Boolean,
                               extras: Option[Extras],
                               parameterLinkSerializer: ParameterLinkSerializer,
                               useManifests: Boolean,
                               complexPathValues: Boolean,
                               instanceTypeSelection: InstanceTypeSelection.InstanceTypeSelection,
                               defaultInstanceType: Option[String],
                               fileResolver: FileSourceResolver,
                               dxApi: DxApi = DxApi.get,
                               logger: Logger = Logger.get,
                               project: DxProject,
                               folder: String)
    extends ExecutableCompiler(extras,
                               parameterLinkSerializer,
                               complexPathValues,
                               fileResolver,
                               dxApi) {
  // database of available instance types for the user/org that owns the project
  private lazy val instanceTypeDb = InstanceType.createDb(Some(project))

  // renderer for job script templates
  private lazy val renderer = Renderer()

  // Preamble required for accessing a private docker registry (if required)
  private lazy val dockerRegistry: Option[DockerRegistry] = extras.flatMap(_.dockerRegistry)
  private lazy val dockerPreamble: Option[String] = {
    def checkFile(uri: String, name: String): Unit = {
      try {
        logger.ignore(dxApi.resolveFile(uri))
      } catch {
        case e: Throwable =>
          throw new Exception(s"""|${name} has to point to a platform file.
                                  |It is now:
                                  |   ${uri}
                                  |Error:
                                  |  ${e}
                                  |""".stripMargin)
      }
    }
    dockerRegistry.map {
      case DockerRegistry(registry, credentials, Some(username), None) =>
        // check that the credentials file is a valid platform path
        checkFile(credentials, ApplicationCompiler.CredentialsKey)
        // render the preamble
        renderer.render(
            ApplicationCompiler.GenericDockerPreambleTemplate,
            Map(
                ApplicationCompiler.RegistryKey -> registry,
                ApplicationCompiler.UsernameKey -> username,
                // strip the URL from the dx:// prefix, so we can use dx-download directly
                ApplicationCompiler.CredentialsKey -> credentials.substring(
                    DxPath.DxUriPrefix.length
                )
            )
        )
      case DockerRegistry(registry, credentials, None, Some(awsRegion)) =>
        checkFile(credentials, ApplicationCompiler.CredentialsKey)
        renderer.render(
            ApplicationCompiler.EcrDockerPreambleTemplate,
            Map(
                ApplicationCompiler.RegistryKey -> registry,
                // strip the URL from the dx:// prefix, so we can use dx-download directly
                ApplicationCompiler.CredentialsKey -> credentials.substring(
                    DxPath.DxUriPrefix.length
                ),
                ApplicationCompiler.AwsRegionKey -> awsRegion
            )
        )
      case _ =>
        throw new Exception(s"invalid dockerRegistry settings ${dockerRegistry}")
    }
  }

  private def generateJobScript(applet: Application): String = {
    val templateAttrs: Map[String, Any] = Map(
        "runtimeJar" -> runtimeJar,
        "runtimeTraceLevel" -> runtimeTraceLevel,
        "streamFiles" -> streamFiles,
        "waitOnUpload" -> waitOnUpload
    )
    applet.kind match {
      case ExecutableKindApplet =>
        val template = applet.instanceType match {
          case DynamicInstanceType => ApplicationCompiler.DynamicAppletCommandTemplate
          case _                   => ApplicationCompiler.StaticAppletCommandTemplate
        }
        renderer.render(
            template,
            templateAttrs ++ Map(
                "dockerPreamble" -> dockerPreamble,
                "runtimePathConfig" -> runtimePathConfig
            )
        )
      case _: ExecutableKindWfFragment =>
        renderer.render(
            ApplicationCompiler.WorkflowFragmentTemplate,
            templateAttrs ++ Map(
                "separateOutputs" -> separateOutputs,
                "waitOnUpload" -> waitOnUpload,
                "runtimePathConfig" -> runtimePathConfig
            )
        )
      case other =>
        ExecutableKind.getCommand(other) match {
          case Some(command) =>
            renderer.render(
                ApplicationCompiler.WorkflowCommandTemplate,
                templateAttrs ++ Map(
                    "command" -> command,
                    "separateOutputs" -> separateOutputs,
                    "waitOnUpload" -> waitOnUpload,
                    "runtimePathConfig" -> runtimePathConfig
                )
            )
          case _ =>
            throw new RuntimeException(
                s"should not generate job script for kind ${other}"
            )
        }
    }
  }

  private def createRunSpec(
      applet: Application,
      executableDict: Map[String, ExecutableLink]
  ): (JsValue, Map[String, JsValue]) = {
    val instanceType: String = applet.instanceType match {
      case static: StaticInstanceType                => instanceTypeDb.apply(static.req).name
      case DefaultInstanceType | DynamicInstanceType =>
        // TODO: should we use the project default here rather than picking one from the database?
        defaultInstanceType.getOrElse(instanceTypeDb.defaultInstanceType.name)
    }
    // Generate the applet's job script
    val jobScript = generateJobScript(applet)
    // build the run spec
    val runSpecRequired = Map(
        "code" -> JsString(jobScript),
        "interpreter" -> JsString("bash"),
        "systemRequirements" ->
          JsObject(
              "main" ->
                JsObject("instanceType" -> JsString(instanceType))
          ),
        "distribution" -> JsString(Constants.OsDistribution),
        "release" -> JsString(Constants.OsRelease),
        "version" -> JsString(Constants.OsVersion)
    )
    // Add default timeout
    val defaultTimeout =
      DxRunSpec.toApiJson(
          DxRunSpec(
              access = None,
              executionPolicy = None,
              restartableEntryPoints = None,
              timeoutPolicy = Some(
                  DxTimeout(Some(ApplicationCompiler.DefaultAppletTimeoutInDays), Some(0), Some(0))
              )
          )
      )
    // Start with the default dx-attribute section, and override
    // any field that is specified in the runtime hints or the individual task section.
    val extrasOverrides =
      extras.flatMap(_.defaultTaskDxAttributes).map(_.getApiRunSpecJson).getOrElse(Map.empty)
    // runtime hints in the task override defaults from extras
    val taskOverrides: Map[String, JsValue] = applet.requirements.collect {
      case RestartRequirement(maxRestarts, default, errors) =>
        val defaultMap: Map[String, Long] = default match {
          case Some(i) => Map("*" -> i)
          case _       => Map.empty
        }
        val restartOn = errors ++ defaultMap match {
          case m if m.isEmpty => None
          case m              => Some(m)
        }
        DxRunSpec.createApiExecutionPolicy(restartOn, maxRestarts)
      case TimeoutRequirement(days, hours, minutes) =>
        DxRunSpec.createApiTimeoutPolicy(days, hours, minutes)
    }.toMap
    // task-specific settings from extras override runtime hints in the task
    val taskSpecificOverrides = applet.kind match {
      case ExecutableKindApplet =>
        extras.flatMap(_.perTaskDxAttributes.flatMap(_.get(applet.name))) match {
          case Some(dta) => dta.getApiRunSpecJson
          case None      => Map.empty
        }
      case _ => Map.empty
    }
    // Add platform Docker image dependency to bundledDepends
    val bundledDependsDocker: Option[JsValue] = applet.container match {
      case DxFileDockerImage(_, dxfile) => {
        val id = dxfile.id

        // If source project not specified for Docker image file, try to find it
        val sourceProject = dxfile.project.getOrElse(
            Try {
              dxApi.getObject(id).describe() match {
                case f: DxFileDescribe => dxApi.project(f.project)
                case _                 => throw new RuntimeException(s"Expected ${id} to be a file")
              }
            } match {
              case Success(project) => project
              case Failure(ex)      => throw new RuntimeException(s"Unable to locate file ${id}", ex)
            }
        )

        // If dependency on Docker image file in another project, clone
        dxApi.cloneDataObject(id, sourceProject, project, folder)

        val mapping = JsObject(
            Constants.BundledDependsNameKey -> JsString(dxfile.describe().name),
            Constants.BundledDependsIdKey -> JsObject(DxUtils.DxLinkKey -> JsString(dxfile.id)),
            Constants.BundledDependsStagesKey -> JsArray(Vector.empty)
        )
        Some(mapping)
      }
      case _ => None
    }

    // Add executables called by this applet to bundledDepends to ensure cloning
    val bundledDependsExec: Vector[JsValue] = executableDict
      .collect {
        case (name, ExecutableLink(_, _, _, applet: DxApplet))     => (name, applet.id)
        case (name, ExecutableLink(_, _, _, workflow: DxWorkflow)) => (name, workflow.id)
      }
      .map {
        case (name, id) =>
          JsObject(
              Constants.BundledDependsNameKey -> JsString(name),
              Constants.BundledDependsIdKey -> JsObject(DxUtils.DxLinkKey -> JsString(id)),
              Constants.BundledDependsStagesKey -> JsArray(Vector.empty)
          )
      }
      .toVector

    // Include runtimeAsset in bundledDepends
    val bundledDepends =
      Vector(runtimeAsset, bundledDependsDocker, bundledDependsExec).flatten match {
        case Vector()            => Map.empty
        case bundledDependsItems => Map(Constants.BundledDependsKey -> JsArray(bundledDependsItems))
      }

    val runSpec = JsObject(
        runSpecRequired ++
          defaultTimeout ++
          extrasOverrides ++
          taskOverrides ++
          taskSpecificOverrides ++
          bundledDepends
    )
    // Add hard-coded instance type info to details
    val instanceTypeDetails: Map[String, JsValue] = applet.instanceType match {
      case StaticInstanceType(
          InstanceTypeRequest(Some(staticInstanceType), _, _, _, _, _, _, _, _, _, _)
          ) =>
        Map(Constants.StaticInstanceType -> JsString(staticInstanceType))
      case _ => Map.empty
    }
    // Add Docker image info to details
    val dockerImageDetails: Map[String, JsValue] = applet.container match {
      case DxFileDockerImage(_, dxfile) => Map(Constants.DockerImage -> dxfile.asJson)
      case NetworkDockerImage(pullName) => Map(Constants.NetworkDockerImage -> JsString(pullName))
      case DynamicDockerImage           => Map(Constants.DynamicDockerImage -> JsBoolean(true))
      case NoImage                      => Map.empty
    }
    (runSpec, instanceTypeDetails ++ dockerImageDetails)
  }

  // Convert the applet meta to JSON, and overlay details from default and task-specific extras
  private def applicationAttributesToNative(
      applet: Application,
      defaultTags: Set[String],
      extendedDescription: Option[String]
  ): (Map[String, JsValue], Map[String, JsValue]) = {
    // Default attributes from extras
    val (defaultMeta, defaultDetails) =
      extras
        .flatMap(_.defaultTaskDxAttributes)
        .map { attr =>
          (attr.getMetaJson, attr.getDetailsJson)
        }
        .getOrElse((Map.empty[String, JsValue], Map.empty[String, JsValue]))
    // Attributes specified in WDL meta section
    val (commonMeta, commonDetails) = callableAttributesToNative(
        callable = applet,
        defaultTags = defaultTags,
        extendedDescription = extendedDescription
    )
    val applicationMeta = applet.attributes.collect {
      case DeveloperNotesAttribute(text)     => "developerNotes" -> JsString(text)
      case VersionAttribute(text)            => "version" -> JsString(text)
      case OpenSourceAttribute(isOpenSource) => "openSource" -> JsBoolean(isOpenSource)
      case CategoriesAttribute(categories)   => "categories" -> JsArray(categories.map(JsString(_)))
    }
    // Defaults and those specified in WDL meta can be overridden by task-specific extras
    val (taskSpecificMeta, taskSpecificDetails) = Option
      .when(applet.kind == ExecutableKindApplet) {
        extras
          .flatMap(
              _.perTaskDxAttributes.flatMap(_.get(applet.name))
          )
          .map { attr =>
            (attr.getMetaJson, attr.getDetailsJson)
          }
      }
      .flatten
      .getOrElse((Map.empty[String, JsValue], Map.empty[String, JsValue]))
    (defaultMeta ++ commonMeta ++ applicationMeta ++ taskSpecificMeta,
     defaultDetails ++ commonDetails ++ taskSpecificDetails)
  }

  // TODO: Use SSP templates for Markdown dependency report

  // Summarize dependencies that are not bundled for cloning
  // with the app/let for adding to description in Markdown
  private def summarizeReportableDependencies(
      fileDependencies: Vector[String],
      runSpecDetails: Map[String, JsValue]
  ): Option[String] = {
    val headerMd =
      "# This app requires access to the following dependencies that are not " +
        "packaged with the app"
    val filesMd = fileDependencies.map(listMd2).mkString match {
      case s: String if s.nonEmpty => s"${listMd("Files")}${s}"
      case _                       => ""
    }
    val networkDockerImageMd = runSpecDetails.get(Constants.NetworkDockerImage) match {
      case Some(JsString(s)) => s"${listMd("Network Docker image")}${listMd2(s)}"
      case _                 => ""
    }
    val dynamicDockerImageMd = runSpecDetails.get(Constants.DynamicDockerImage) match {
      case Some(JsBoolean(true)) => listMd("Docker image determined at runtime")
      case _                     => ""
    }
    val staticInstanceTypeMd = runSpecDetails.get(Constants.StaticInstanceType) match {
      case Some(JsString(s)) => s"${listMd("Hard-coded instance type")}${listMd2(s)}"
      case _                 => ""
    }
    val staticInstanceTypeSelectionMd = instanceTypeSelection match {
      case InstanceTypeSelection.Static =>
        listMd(
            "Available instance types determined during WDL compilation will be " +
              "used to select instance types during workflow execution"
        )
      case _ => ""
    }
    val md =
      s"${filesMd}${networkDockerImageMd}${dynamicDockerImageMd}${staticInstanceTypeMd}" +
        s"${staticInstanceTypeSelectionMd}"
    Option.when(md.nonEmpty)(s"${headerMd}${md}")
  }

  def createAccess(applet: Application): JsValue = {
    // defaults are taken from
    // extras global defaults < task runtime section < task-specific extras
    val defaultAccess: DxAccess = extras.map(_.getDefaultAccess).getOrElse(DxAccess.empty)
    val taskAccess: DxAccess = applet.requirements
      .collectFirst {
        case AccessRequirement(network, project, allProjects, developer, projectCreation) =>
          val networkOpt = if (network.isEmpty) {
            None
          } else {
            Some(network)
          }
          DxAccess(networkOpt,
                   project.map(DxAccessLevel.withName),
                   allProjects.map(DxAccessLevel.withName),
                   developer,
                   projectCreation)
      }
      .getOrElse(DxAccess.empty)
    val taskSpecificAccess: DxAccess = (applet.kind match {
      case ExecutableKindApplet =>
        extras.map(_.getTaskAccess(applet.name))
      case _ => None
    }).getOrElse(DxAccess.empty)
    // If we are using a private docker registry, add the allProjects: VIEW
    // access to tasks.
    val allProjectsAccess: DxAccess = dockerRegistry match {
      case None    => DxAccess.empty
      case Some(_) => DxAccess.empty.copy(allProjects = Some(DxAccessLevel.View))
    }
    // update depending on applet type
    val appletKindAccess = applet.kind match {
      case ExecutableKindApplet =>
        applet.container match {
          // Require network access if Docker image needs to be downloaded
          case NetworkDockerImage(_) | DynamicDockerImage =>
            Some(DxAccess.empty.copy(network = Some(Vector("*"))))
          case _ => None
        }
      case ExecutableKindWorkflowOutputReorg =>
        // The reorg applet requires higher permissions to organize the output directory.
        Some(DxAccess.empty.copy(project = Some(DxAccessLevel.Contribute)))
      case _ =>
        // Scatters need network access, because they spawn subjobs that (may) use dx-docker.
        // We end up allowing all applets to use the network
        Some(DxAccess.empty.copy(network = Some(Vector("*"))))
    }
    // using manifests requires at least UPLOAD access
    val manifestAccess = if (useManifests) {
      Some(DxAccess.empty.copy(project = Some(DxAccessLevel.Upload)))
    } else {
      None
    }
    // merge all
    val access = defaultAccess
      .merge(taskAccess)
      .merge(taskSpecificAccess)
      .merge(allProjectsAccess)
      .mergeOpt(appletKindAccess)
      .mergeOpt(manifestAccess)
    access.toJson match {
      case JsObject(fields) if fields.isEmpty => JsNull
      case fields                             => fields
    }
  }

  /**
    * Builds an '/applet/new' request.
    * For applets that call other applets, we pass a directory of the callees,
    * so they can be found at runtime.
    * @param applet applet IR
    * @param executableDict mapping of callable names to executables
    * @return
    */
  def apply(
      applet: Application,
      executableDict: Map[String, ExecutableLink]
  ): Map[String, JsValue] = {
    logger.trace(s"Building /applet/new request for ${applet.name}")
    // convert inputs and outputs to dxapp inputSpec, automatically add requirements/hints parameters
    val inputParams = if (useManifests) {
      Vector(
          ExecutableCompiler.InputManifestParameter,
          ExecutableCompiler.InputManfestFilesParameter,
          ExecutableCompiler.InputLinksParameter,
          ExecutableCompiler.WorkflowInputManifestParameter,
          ExecutableCompiler.WorkflowInputManfestFilesParameter,
          ExecutableCompiler.WorkflowInputLinksParameter,
          ExecutableCompiler.OutputIdParameter,
          ExecutableCompiler.ExtraOutputsParameter,
          ExecutableCompiler.CallNameParameter,
          ExecutableCompiler.OverridesParameter
      )
    } else {
      applet.inputs ++ Vector(
          ExecutableCompiler.OverridesParameter
      )
    }
    val inputSpec = inputParams
      .sortWith(_.name < _.name)
      .flatMap { param =>
        try {
          inputParameterToNative(param)
        } catch {
          case ex: Throwable =>
            throw new Exception(
                s"Error converting input parameter ${param} to native type",
                ex
            )
        }
      }
    // Collect file dependencies from inputs & private variables
    val fileDependencies = (
        applet.inputs
          .filter(param => param.dxType == Type.TFile)
          .flatMap(fileDependenciesFromParam) ++ applet.staticFileDependencies
    ).distinct
    val fileDependenciesDetails = Option
      .when(fileDependencies.nonEmpty)(
          Map(Constants.FileDependencies -> JsArray(fileDependencies.map(s => JsString(s))))
      )
      .getOrElse(Map.empty)
    // convert outputs to outputSpec
    val outputParams = if (useManifests) {
      Vector(Parameter(Constants.OutputManifest, Type.TFile))
    } else {
      applet.outputs
    }
    val outputSpec: Vector[JsValue] = outputParams
      .sortWith(_.name < _.name)
      .flatMap { param =>
        try {
          outputParameterToNative(param)
        } catch {
          case ex: Throwable =>
            throw new Exception(
                s"Error converting output parameter ${param} to native type",
                ex
            )
        }
      }
    // build the dxapp runSpec
    val (runSpec, runSpecDetails) = createRunSpec(applet, executableDict)
    // A fragment is hidden, not visible under default settings. This
    // allows the workflow copying code to traverse it, and link to
    // anything it calls.
    val hidden: Boolean =
      applet.kind match {
        case _: ExecutableKindWfFragment => true
        case _                           => false
      }

    // Create linking information: results in two maps, one that's added to the
    // application details, and one with links to applets that could get called
    // at runtime (if this applet is copied, we need to maintain referential integrity)
    // Note: this doesn't seem to ensure cloning when copying workflow.
    val (dxLinks, linkInfo) = executableDict.map {
      case (name, link) =>
        val linkName = s"link_${name}"
        (
            linkName -> JsObject(DxUtils.DxLinkKey -> JsString(link.dxExec.id)),
            name -> ExecutableLink.serialize(link)
        )
    }.unzip
    // build the details JSON
    val defaultTags = Set(Constants.CompilerTag)
    // Build extended description with dependency report
    val extendedDescription = summarizeReportableDependencies(
        fileDependencies = fileDependencies,
        runSpecDetails = runSpecDetails
    )
    val (taskMeta, taskDetails) =
      applicationAttributesToNative(applet, defaultTags, extendedDescription)
    val delayDetails = delayWorkspaceDestructionToNative
    // meta information used for running workflow fragments
    val metaDetails: Map[String, JsValue] =
      applet.kind match {
        case ExecutableKindWfFragment(_, blockPath, inputs, scatterChunkSize) =>
          Map(
              Constants.ExecLinkInfo -> JsObject(linkInfo.toMap),
              Constants.BlockPath -> JsArray(blockPath.map(JsNumber(_))),
              Constants.WfFragmentInputTypes -> TypeSerde.serializeSpec(inputs.map {
                case (dxName, t) => dxName.decoded -> t
              })
          ) ++ scatterChunkSize
            .map(chunkSize => Map(Constants.ScatterChunkSize -> JsNumber(chunkSize)))
            .getOrElse(Map.empty)
        case ExecutableKindWfInputs(blockPath) if blockPath.nonEmpty =>
          val types = applet.inputVars.map(p => p.name -> p.dxType).toMap
          Map(
              Constants.BlockPath -> JsArray(blockPath.map(JsNumber(_))),
              Constants.WfFragmentInputTypes -> TypeSerde.serializeSpec(types.map {
                case (dxName, t) => dxName.decoded -> t
              })
          )
        case ExecutableKindWfOutputs(blockPath) if blockPath.nonEmpty =>
          val types = applet.inputVars.map(p => p.name -> p.dxType).toMap
          Map(
              Constants.BlockPath -> JsArray(blockPath.map(JsNumber(_))),
              Constants.WfFragmentInputTypes -> TypeSerde.serializeSpec(types.map {
                case (dxName, t) => dxName.decoded -> t
              })
          )
        case _: ExecutableKindWfInputs | _: ExecutableKindWfOutputs |
            ExecutableKindWfCustomReorgOutputs | ExecutableKindWorkflowOutputReorg =>
          val types = applet.inputVars.map(p => p.name -> p.dxType).toMap
          Map(Constants.WfFragmentInputTypes -> TypeSerde.serializeSpec(types.map {
            case (dxName, t) => dxName.decoded -> t
          }))
        case _ => Map.empty
      }
    // compress and base64 encode the source code
    val sourceDocStringEncoded = CodecUtils.gzipAndBase64Encode(applet.document.getDocContents)
    // compress and base64 encode the instance types, unless we specify that we want to
    // resolve them at runtime, which requires that the user running the applet has
    // permission to describe the project it is running in
    val dbEncoded = Option.when(instanceTypeSelection == InstanceTypeSelection.Static) {
      CodecUtils.gzipAndBase64Encode(instanceTypeDb.toJson.prettyPrint)
    }
    // serilize default runtime attributes
    val defaultRuntimeAttributes = extras
      .flatMap(ex =>
        ex.defaultRuntimeAttributes.map(attr => JsObject(ValueSerde.serializeMap(attr)))
      )
    val auxDetails = Vector(
        Some(Constants.DocContents -> JsString(sourceDocStringEncoded)),
        Some(Constants.ParseOptions -> applet.document.optionsToJson),
        dbEncoded.map(db => Constants.InstanceTypeDb -> JsString(db)),
        defaultRuntimeAttributes.map(attr => Constants.RuntimeAttributes -> attr),
        Option.when(useManifests)(Constants.UseManifests -> JsBoolean(true)),
        Option.when(complexPathValues)(Constants.PathsAsObjects -> JsBoolean(true))
    ).flatten.toMap
    // combine all details into a single Map
    val details: Map[String, JsValue] =
      taskDetails ++
        runSpecDetails ++
        delayDetails ++
        dxLinks.toMap ++
        metaDetails ++
        auxDetails ++
        fileDependenciesDetails
    // build the API request
    val requestRequired = Map(
        "name" -> JsString(applet.dxName),
        "inputSpec" -> JsArray(inputSpec),
        "outputSpec" -> JsArray(outputSpec),
        "runSpec" -> runSpec,
        "dxapi" -> JsString(dxApi.version),
        "details" -> JsObject(details),
        "hidden" -> JsBoolean(hidden)
    )
    // look for ignoreReuse in runtime hints and in extras - the later overrides the former
    val ignoreReuse = applet.requirements
      .collectFirst {
        case IgnoreReuseRequirement(value) => value
      }
      .orElse(
          extras.flatMap(_.ignoreReuse)
      )
      .map(ignoreReuse => Map("ignoreReuse" -> JsBoolean(ignoreReuse)))
      .getOrElse(Map.empty)
    // build the dxapp access section
    val access = createAccess(applet) match {
      case JsNull  => Map.empty
      case jsValue => Map("access" -> jsValue)
    }
    taskMeta ++ requestRequired ++ access ++ ignoreReuse
  }
}
