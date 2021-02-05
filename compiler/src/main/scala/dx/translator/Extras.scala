package dx.translator

// Place to put any extra options, equivalent to Cromwell workflowOptions.
// Also, allows dnanexus specific configuration per task.

import java.nio.file.Path
import dx.api._
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
    val Access = "access"
    val ExecutionPolicy = "executionPolicy"
    val RestartableEntryPoints = "restartableEntryPoints"
    val TimeoutPolicy = "timeoutPolicy"
    val EntryPointNames = Set("all", "master")

    def read(jsv: JsValue): DxRunSpec = {
      jsv match {
        case JsObject(fields) =>
          val restartableEntryPoints =
            JsUtils.getOptionalString(fields, RestartableEntryPoints).map {
              case name if EntryPointNames.contains(name) => name
              case name =>
                throw new Exception(s"Unsupported restartableEntryPoints value ${name}")
            }
          val timeout = JsUtils.getOptionalFields(fields, TimeoutPolicy).map { timeout =>
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
          DxRunSpec(
              fields.get(Access).map(_.convertTo[DxAccess]),
              fields.get(ExecutionPolicy).map(_.convertTo[DxExecPolicy]),
              restartableEntryPoints,
              timeout
          )
        case _ =>
          deserializationError(s"invalid runSpec ${jsv}")
      }
    }

    def write(runSpec: DxRunSpec): JsValue = {
      val fields = Vector(
          runSpec.access.map(x => Access -> x.toJson),
          runSpec.timeoutPolicy.map(x => TimeoutPolicy -> JsObject("*" -> x.toJson)),
          runSpec.execPolicy.map(x => ExecutionPolicy -> x.toJson),
          runSpec.restartableEntryPoints.map(x => RestartableEntryPoints -> JsString(x))
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
        case other =>
          deserializationError(s"invalid or missing appUri value ${other}")
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
        case Some(JsString(uri)) if uri.trim.startsWith(DxPath.DxUriPrefix) =>
          // if provided, check that the fileID is valid and present
          // format dx file ID
          val reorgFileID: String = uri.trim.replace(DxPath.DxUriPrefix, "")
          // if input file ID is invalid, DxFile.getInstance will thow an IllegalArgumentException
          // if reorgFileID cannot be found, describe will throw a ResourceNotFoundException
          Logger.get.ignore(DxApi.get.file(reorgFileID).describe())
          Some(uri)
        case Some(JsString(uri)) if uri.trim.isEmpty => None
        case _ =>
          deserializationError(
              """In the 'custom_reorg' section of extras, 'configFile' must be specified as
                |a valid DNAnexus file in the form 'dx://file-XXX'. Please set the value
                |to null if there is no configuration file.""".stripMargin.replaceAll("\n", "")
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
  implicit val execPolicyFormat: RootJsonFormat[DxExecPolicy] = jsonFormat2(DxExecPolicy)
  implicit val timeoutPolicyFormat: RootJsonFormat[DxTimeout] = jsonFormat3(DxTimeout)
  implicit val licenseFormat: RootJsonFormat[DxLicense] = jsonFormat6(DxLicense)
  implicit val detailsFormat: RootJsonFormat[DxDetails] = jsonFormat1(DxDetails)
  implicit val dxAppFormat: RootJsonFormat[DxAppJson] = jsonFormat2(DxAppJson)
  implicit val scatterAttrsFormat: RootJsonFormat[DxScatterAttrs] = jsonFormat1(DxScatterAttrs)
  implicit val workflowAttrsFormat: RootJsonFormat[DxWorkflowAttrs] = jsonFormat2(DxWorkflowAttrs)
  implicit val dockerRegistryFormat: RootJsonFormat[DockerRegistry] = jsonFormat5(DockerRegistry)
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

case class DxExecPolicy(restartOn: Option[Map[String, Long]], maxRestarts: Option[Long])

case class DxTimeout(days: Option[Long], hours: Option[Long], minutes: Option[Long])

case class DxRunSpec(access: Option[DxAccess],
                     execPolicy: Option[DxExecPolicy],
                     restartableEntryPoints: Option[String],
                     timeoutPolicy: Option[DxTimeout])

object DxRunSpec {
  def toApiJson(runSpec: DxRunSpec): Map[String, JsValue] = {
    runSpec.toJson.asJsObject.fields.filterNot {
      // the access field is in runSpec in in dxapp.json but not in the API call
      case (key, _) => key == "access"
    }
  }
}

case class DxLicense(name: String,
                     repoUrl: String,
                     version: String,
                     license: String,
                     licenseUrl: String,
                     author: String)

case class DxDetails(upstreamProjects: Option[Vector[DxLicense]])

case class DxAppJson(runSpec: Option[DxRunSpec], details: Option[DxDetails]) {
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
case class DxScatterAttrs(chunkSize: Option[Int] = None)

/**
  * Runtime attributes set at a per-workflow level.
  * @param scatterDefaults default scatter attributes that apply to an entire workflow
  * @param perScatterAttrs scatter attributes that apply to individual scatter blocks
  */
case class DxWorkflowAttrs(scatterDefaults: Option[DxScatterAttrs],
                           perScatterAttrs: Option[Map[String, DxScatterAttrs]])

case class DockerRegistry(registry: String,
                          credentials: String,
                          username: Option[String],
                          awsConfig: Option[String],
                          awsProfile: Option[String])

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
                  perWorkflowDxAttributes: Option[Map[String, DxWorkflowAttrs]],
                  dockerRegistry: Option[DockerRegistry],
                  customReorgAttributes: Option[CustomReorgSettings],
                  ignoreReuse: Option[Boolean],
                  delayWorkspaceDestruction: Option[Boolean]) {
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
    // the format use to have some snake-cased keys, so we update them all to camel-case first
    camelizeKeys(jsv).convertTo[Extras]
  }

  def parse(path: Path): Extras = {
    parse(JsUtils.jsFromFile(path))
  }
}
