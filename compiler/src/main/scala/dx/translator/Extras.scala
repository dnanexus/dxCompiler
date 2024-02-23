package dx.translator

// Place to put any extra options, equivalent to Cromwell workflowOptions.
// Also, allows DNAnexus specific configuration per task.

import java.nio.file.Path
import dx.api._
import dx.core.Constants
import dx.core.ir.{Value, ValueSerde}
import dx.util.{JsUtils, Logger}
import spray.json._

object ExtrasJsonProtocol extends DefaultJsonProtocol {
  implicit object DxAccessLevelProtocol extends RootJsonFormat[DxAccessLevel] {
    def write(level: DxAccessLevel): JsString = JsString(level.name)
    def read(jsv: JsValue): DxAccessLevel = {
      jsv match {
        case JsString(name) => DxAccessLevel.withName(name)
        case other          => deserializationError(s"invalid access level ${other}")
      }
    }
  }

  implicit object DxRunSpecJsonProtocol extends RootJsonFormat[DxRunSpec] {
    def read(jsv: JsValue): DxRunSpec = {
      jsv match {
        case JsObject(fields) =>
          val restartableEntryPoints =
            JsUtils.getOptionalString(fields, DxRunSpec.RestartableEntryPoints).map {
              case name if DxRunSpec.EntryPointNames.contains(name) => name
              case name =>
                throw new Exception(s"Unsupported restartableEntryPoints value ${name}")
            }
          val timeout = JsUtils.getOptionalFields(fields, DxRunSpec.TimeoutPolicy).map { timeout =>
            if (timeout.size != 1) {
              deserializationError("Exactly one entry-point timeout can be specified")
            }
            if (timeout.keys.head != "*") {
              deserializationError(
                  """Only a general timeout for all entry points is supported ("*")"""
              )
            }
            timeout("*").convertTo[DxTimeout]
          }
          val headJobOnDemand = JsUtils.getOptionalBoolean(fields, DxRunSpec.HeadJobOnDemand)
          DxRunSpec(
              fields.get(DxRunSpec.Access).map(_.convertTo[DxAccess]),
              fields.get(DxRunSpec.ExecutionPolicy).map(_.convertTo[DxExecPolicy]),
              restartableEntryPoints,
              timeout,
              headJobOnDemand
          )
        case _ =>
          deserializationError(s"invalid runSpec ${jsv}")
      }
    }

    def write(runSpec: DxRunSpec): JsValue = {
      val fields = Vector(
          runSpec.access.map(x => DxRunSpec.Access -> x.toJson),
          runSpec.timeoutPolicy.map(x => DxRunSpec.TimeoutPolicy -> JsObject("*" -> x.toJson)),
          runSpec.executionPolicy.map(x => DxRunSpec.ExecutionPolicy -> x.toJson),
          runSpec.restartableEntryPoints.map(x => DxRunSpec.RestartableEntryPoints -> JsString(x)),
          runSpec.headJobOnDemand.map(x => DxRunSpec.HeadJobOnDemand -> JsBoolean(x))
      ).flatten.toMap
      if (fields.isEmpty) {
        JsNull
      } else {
        JsObject(fields)
      }
    }
  }

  implicit object IrValueFormat extends RootJsonFormat[Value] {
    def read(jsv: JsValue): Value = {
      ValueSerde.deserialize(jsv)
    }

    def write(value: Value): JsValue = {
      ValueSerde.serialize(value)
    }
  }

  implicit object CustomReorgFormat extends RootJsonFormat[CustomReorgSettings] {
    val AppUri = "appUri"
    val ConfigFile = "configFile"

    def read(jsv: JsValue): CustomReorgSettings = {
      val fields = jsv.asJsObject.fields
      val reorgAppId: String = fields.get(AppUri).orElse(fields.get("app_id")) match {
        case Some(JsString(uri)) => uri
        case None =>
          deserializationError(s"appUri must be specified in the customReorgAttributes section")
        case other =>
          deserializationError(s"invalid appUri value ${other}")
      }
      val access: Option[JsValue] =
        try {
          DxApi.get.getObject(reorgAppId) match {
            case exe: DxExecutable =>
              exe.describe(Set(Field.Access)) match {
                case desc: DxAppDescribe    => desc.access
                case desc: DxAppletDescribe => desc.access
                case _ =>
                  deserializationError(s"${reorgAppId} is not a DNAnexus executable")
              }
          }
        } catch {
          case t: Throwable =>
            deserializationError(s"Invalid reorg app(let) ID ${reorgAppId}", t)
        }
      val hasAccess = access match {
        case Some(JsObject(x)) =>
          JsObject(x).fields.get("project") match {
            case Some(JsString("CONTRIBUTE")) | Some(JsString("ADMINISTER")) => true
            case _                                                           => false
          }
        case _ => false
      }
      if (!hasAccess) {
        throw new dx.PermissionDeniedException(
            s"""ERROR: App(let) for custom reorg stage ${reorgAppId} does not
               |have CONTRIBUTOR or ADMINISTRATOR access and this is required.""".stripMargin
              .replaceAll("\n", "")
        )
      }
      val reorgConf: Option[String] = fields.get(ConfigFile).orElse(fields.get("conf")) match {
        case Some(JsString(uri)) if uri.trim.isEmpty                        => None
        case Some(JsString(uri)) if uri.trim.startsWith(DxPath.DxUriPrefix) =>
          // if provided, check that the fileID is valid and present format dx file ID
          val dxfile = DxApi.get.resolveFile(uri.trim)
          // if input file ID is invalid, DxFile.getInstance will thow an IllegalArgumentException
          // if reorgFileID cannot be found, describe will throw a ResourceNotFoundException
          Logger.get.ignore(dxfile.describe())
          Some(dxfile.asUri)
        case Some(JsNull) => None
        case _ =>
          deserializationError(
              """In the 'customReorgAttributes' section of extras, 'configFile' must be specified as
                |a valid DNAnexus file in the form 'dx://file-XXX'. Please set the value
                |to null if there is no configuration file.""".stripMargin.replaceAll("\n", " ")
          )
      }
      Logger.get.trace(
          s"""|Writing your own applet for reorganization purposes is tricky. If you are not careful,
              |it may misplace or outright delete files.
              |
              |The applet ${reorgAppId} requires CONTRIBUTE project access so it can move files and
              |folders. The applet must be idempotent, so that if the instance it runs on crashes,
              |it can safely restart. It must also be careful about inputs that are also outputs:
              |generally, these should not be moved. It should use bulk object operations, so as
              |not to overload the API server.
              |
              |You may refer to this example:
              |https://github.com/dnanexus/dxCompiler/blob/master/doc/ExpertOptions.md#use-your-own-applet
            """.stripMargin.replaceAll("\n", " ")
      )
      CustomReorgSettings(reorgAppId, reorgConf)
    }

    override def write(obj: CustomReorgSettings): JsValue = {
      JsObject(
          Vector(
              Some(AppUri -> JsString(obj.appUri)),
              obj.configFile.map(uri => ConfigFile -> JsString(uri))
          ).flatten.toMap
      )
    }
  }

  implicit val accessFormat: RootJsonFormat[DxAccess] = jsonFormat5(DxAccess.apply)
  implicit val executionPolicyFormat: RootJsonFormat[DxExecPolicy] = jsonFormat2(DxExecPolicy.apply)
  implicit val timeoutPolicyFormat: RootJsonFormat[DxTimeout] = jsonFormat3(DxTimeout)
  implicit val licenseFormat: RootJsonFormat[DxLicense] = jsonFormat6(DxLicense)
  implicit val detailsFormat: RootJsonFormat[DxDetails] = jsonFormat1(DxDetails)
  implicit val dxAppFormat: RootJsonFormat[DxAppJson] = jsonFormat13(DxAppJson)
  implicit val scatterAttrsFormat: RootJsonFormat[DxScatterAttrs] = jsonFormat1(DxScatterAttrs)
  implicit val workflowAttrsFormat: RootJsonFormat[DxWorkflowAttrs] = jsonFormat12(DxWorkflowAttrs)
  implicit val dockerRegistryFormat: RootJsonFormat[DockerRegistry] = jsonFormat4(DockerRegistry)
  implicit val defaultReorgSettingsFormat: RootJsonFormat[DefaultReorgSettings] = jsonFormat1(
      DefaultReorgSettings
  )
  implicit val extrasFormat: RootJsonFormat[Extras] = jsonFormat8(Extras.apply)
}

import ExtrasJsonProtocol._

object DxAccess {
  val empty: DxAccess = DxAccess(None, None, None, None, None)
}

case class DxAccess(network: Option[Vector[String]],
                    project: Option[DxAccessLevel],
                    allProjects: Option[DxAccessLevel],
                    developer: Option[Boolean],
                    projectCreation: Option[Boolean]) {

  // Merge with an additional set of access requirments.
  // Take the maximum for each field, and merge the network.
  def merge(from: DxAccess): DxAccess = {
    def mergeAccessLevels(al1: Option[DxAccessLevel],
                          al2: Option[DxAccessLevel]): Option[DxAccessLevel] = {
      Vector(al1, al2).flatten match {
        case Vector(x)    => Some(x)
        case Vector(x, y) => Some(Vector(x, y).max)
        case _            => None
      }
    }
    def mergeBooleans(a: Option[Boolean], b: Option[Boolean]): Option[Boolean] = {
      Vector(a, b).flatten match {
        case Vector(x)    => Some(x)
        case Vector(x, y) => Some(x || y)
        case _            => None
      }
    }
    val networkMerged =
      (network.getOrElse(Vector.empty) ++ from.network.getOrElse(Vector.empty)).toSet match {
        // "*" includes all other addresses
        case s if s.isEmpty       => None
        case s if s.contains("*") => Some(Vector("*"))
        case s                    => Some(s.toVector)
      }
    DxAccess(
        networkMerged,
        mergeAccessLevels(project, from.project),
        mergeAccessLevels(allProjects, from.allProjects),
        mergeBooleans(developer, from.developer),
        mergeBooleans(projectCreation, from.projectCreation)
    )
  }

  def mergeOpt(from: Option[DxAccess]): DxAccess = {
    from match {
      case Some(access) => merge(access)
      case None         => this
    }
  }
}

case class DxExecPolicy(restartOn: Option[Map[String, Long]], maxRestarts: Option[Long]) {
  restartOn.foreach { r =>
    val unsupported = r.keySet.diff(DxExecPolicy.RunSpecExecPolicyRestartOnAttrs)
    if (unsupported.nonEmpty) {
      deserializationError(s"unsupported field(s) ${unsupported.mkString(",")} in restart policy")
    }
  }
}

object DxExecPolicy {
  private val RunSpecExecPolicyRestartOnAttrs = Set(
      "ExecutionError",
      "UnresponsiveWorker",
      "JMInternalError",
      "AppInternalError",
      "JobTimeoutExceeded",
      "*"
  )
}

case class DxTimeout(days: Option[Long], hours: Option[Long], minutes: Option[Long])

case class DxRunSpec(access: Option[DxAccess],
                     executionPolicy: Option[DxExecPolicy],
                     restartableEntryPoints: Option[String],
                     timeoutPolicy: Option[DxTimeout],
                     headJobOnDemand: Option[Boolean] = None) {}

object DxRunSpec {
  val Access = "access"
  val ExecutionPolicy = "executionPolicy"
  val RestartableEntryPoints = "restartableEntryPoints"
  val TimeoutPolicy = "timeoutPolicy"
  val EntryPointNames = Set("all", "master")
  val HeadJobOnDemand = "headJobOnDemand"

  def toApiJson(runSpec: DxRunSpec): Map[String, JsValue] = {
    runSpec.toJson.asJsObject.fields.filterNot {
      // the access field is in runSpec in in dxapp.json but not in the API call
      case (key, _) => key == "access"
    }
  }

  def createApiExecutionPolicy(restartOn: Option[Map[String, Long]],
                               maxRestarts: Option[Long]): (String, JsValue) = {
    DxRunSpec.ExecutionPolicy -> DxExecPolicy(restartOn, maxRestarts).toJson
  }

  def createApiTimeoutPolicy(days: Option[Long],
                             hours: Option[Long],
                             minutes: Option[Long]): (String, JsValue) = {
    DxRunSpec.TimeoutPolicy -> JsObject(
        "*" -> DxTimeout(days.orElse(Some(0)), hours.orElse(Some(0)), minutes.orElse(Some(0))).toJson
    )
  }
}

case class DxLicense(name: String,
                     repoUrl: String,
                     version: String,
                     license: String,
                     licenseUrl: String,
                     author: String)

case class DxDetails(upstreamProjects: Option[Vector[DxLicense]])

abstract class DxMeta(title: Option[String],
                      summary: Option[String],
                      description: Option[String],
                      developerNotes: Option[String],
                      version: Option[String],
                      categories: Option[Vector[String]],
                      types: Option[Vector[String]],
                      tags: Option[Vector[String]],
                      properties: Option[Map[String, String]],
                      treeTurnaroundTimeThreshold: Option[Long]) {
  def getMetaJson: Map[String, JsValue] = {
    Vector(
        title.map(t => "title" -> JsString(t)),
        description.map(d => "description" -> JsString(d)),
        summary.map(s => "summary" -> JsString(s)),
        developerNotes.map(d => "developerNotes" -> JsString(d)),
        version.map(v => "version" -> JsString(v)),
        categories.map(c => "categories" -> JsArray(c.map(JsString(_)))),
        types.map(t => "types" -> JsArray(t.map(JsString(_)))),
        tags.map(t => "tags" -> JsArray(t.map(JsString(_)))),
        treeTurnaroundTimeThreshold.map(tat => "treeTurnaroundTimeThreshold" -> JsNumber(tat)),
        properties.map(p =>
          "properties" -> JsObject(p.map {
            case (key, value) => key -> JsString(value)
          })
        )
    ).flatten.toMap
  }
}

case class DxAppJson(runSpec: Option[DxRunSpec] = None,
                     details: Option[DxDetails] = None,
                     title: Option[String] = None,
                     summary: Option[String] = None,
                     description: Option[String] = None,
                     developerNotes: Option[String] = None,
                     version: Option[String] = None,
                     categories: Option[Vector[String]] = None,
                     types: Option[Vector[String]] = None,
                     tags: Option[Vector[String]] = None,
                     properties: Option[Map[String, String]] = None,
                     openSource: Option[Boolean] = None,
                     treeTurnaroundTimeThreshold: Option[Long] = None)
    extends DxMeta(title,
                   summary,
                   description,
                   developerNotes,
                   version,
                   categories,
                   types,
                   tags,
                   properties,
                   treeTurnaroundTimeThreshold) {

  override def getMetaJson: Map[String, JsValue] = {
    super.getMetaJson ++ Vector(
        openSource.map(o => "openSource" -> JsBoolean(o))
    ).flatten.toMap
  }

  def getApiRunSpecJson: Map[String, JsValue] = {
    runSpec.map(DxRunSpec.toApiJson).getOrElse(Map.empty)
  }

  def getDetailsJson: Map[String, JsValue] = {
    details match {
      case None          => Map.empty
      case Some(details) => details.toJson.asJsObject.fields
    }
  }
}

/**
  * Runtime attributes for scatter blocks.
  * @param chunkSize maximum number of scatter jobs to run at once
  */
case class DxScatterAttrs(chunkSize: Option[Int] = None) {
  chunkSize.foreach { size =>
    if (size > Constants.JobsPerScatterLimit) {
      deserializationError(
          s"The number of jobs per scatter must be between 1-${Constants.JobsPerScatterLimit}"
      )
    }
  }
}

/**
  * Runtime attributes set at a per-workflow level.
  * @param scatterDefaults default scatter attributes that apply to an entire workflow
  * @param scatters scatter attributes that apply to individual scatter blocks
  */
case class DxWorkflowAttrs(scatterDefaults: Option[DxScatterAttrs],
                           scatters: Option[Map[String, DxScatterAttrs]],
                           title: Option[String],
                           summary: Option[String],
                           description: Option[String],
                           developerNotes: Option[String],
                           version: Option[String],
                           categories: Option[Vector[String]],
                           types: Option[Vector[String]],
                           tags: Option[Vector[String]],
                           properties: Option[Map[String, String]],
                           treeTurnaroundTimeThreshold: Option[Long])
    extends DxMeta(title,
                   summary,
                   description,
                   developerNotes,
                   version,
                   categories,
                   types,
                   tags,
                   properties,
                   treeTurnaroundTimeThreshold)

case class DockerRegistry(registry: String,
                          credentials: String,
                          username: Option[String],
                          awsRegion: Option[String]) {
  if (!(username.isDefined || awsRegion.isDefined)) {
    deserializationError(
        "either 'username' or 'awsRegion' must be defined in Extras.dockerRegistry"
    )
  }
}

sealed trait ReorgSettings {
  val enabled: Boolean
}
case class DefaultReorgSettings(enabled: Boolean) extends ReorgSettings
case class CustomReorgSettings(appUri: String,
                               configFile: Option[String] = None,
                               enabled: Boolean = true)
    extends ReorgSettings

case class Extras(defaultRuntimeAttributes: Option[Map[String, Value]],
                  defaultTaskDxAttributes: Option[DxAppJson],
                  perTaskDxAttributes: Option[Map[String, DxAppJson]],
                  defaultWorkflowDxAttributes: Option[DxWorkflowAttrs],
                  perWorkflowDxAttributes: Option[Map[String, DxWorkflowAttrs]],
                  dockerRegistry: Option[DockerRegistry],
                  customReorgAttributes: Option[CustomReorgSettings],
                  ignoreReuse: Option[Boolean]) {
  defaultRuntimeAttributes.foreach { attrs =>
    val unsupportedRuntimeAttrs = attrs.keySet.diff(Extras.RuntimeAttrs)
    if (unsupportedRuntimeAttrs.nonEmpty) {
      deserializationError(
          s"""|Unsupported runtime attribute(s) ${unsupportedRuntimeAttrs.mkString(",")};
              |we currently support ${Extras.RuntimeAttrs}
              |""".stripMargin.replaceAll("\n", " ")
      )
    }
  }

  def getDefaultAccess: DxAccess = {
    defaultTaskDxAttributes.flatMap(_.runSpec.flatMap(_.access)).getOrElse(DxAccess.empty)
  }

  def getTaskAccess(taskName: String): DxAccess = {
    perTaskDxAttributes
      .flatMap(_.get(taskName).flatMap(_.runSpec.flatMap(_.access)))
      .getOrElse(DxAccess.empty)
  }
}

object Extras {
  private val RuntimeAttrs =
    Set(
        // DNAnexus-specific runtime attributes
        "dx_instance_type",
        "dx_timeout",
        "dx_restart",
        "dx_ignore_reuse",
        // WDL native runtime attributes
        "memory",
        "disks",
        "cpu",
        "docker",
        "container",
        "gpu",
        "maxRetries",
        "returnCodes",
        // CWL runtime attributes used in ResourceRequirement
        "coresMin",
        "coresMax",
        "ramMin",
        "ramMax",
        "tmpdirMin",
        "tmpdirMax",
        "outdirMin",
        "outdirMax"
    )
  private val camelizeRegexp = "_([a-z\\d])".r

  private def camelize(s: String): String = {
    camelizeRegexp.replaceAllIn(s, { m =>
      m.group(1).toUpperCase()
    })
  }

  private def camelizeKeys(jsv: JsValue): JsValue = {
    jsv match {
      case JsObject(fields) =>
        JsObject(fields.map {
          case (key, value) => camelize(key) -> value
        })
      case JsArray(values) => JsArray(values.map(camelizeKeys))
      case _               => jsv
    }
  }

  def parse(jsv: JsValue): Extras = {
    // The format used to have some snake-cased and some weirdly named attributes,
    // so we update them all to camel-case and fix the names first.
    val fixed = JsObject(jsv.asJsObject.fields.map {
      case (key, value) if Set("custom-reorg", "custom_reorg").contains(key) =>
        "customReorgAttributes" -> value
      case other => other
    })
    val camelized = camelizeKeys(fixed)
    camelized.convertTo[Extras]
  }

  def parse(path: Path): Extras = {
    parse(JsUtils.jsFromFile(path))
  }
}
