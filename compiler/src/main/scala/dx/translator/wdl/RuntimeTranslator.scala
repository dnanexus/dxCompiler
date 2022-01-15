package dx.translator.wdl

import dx.api.{DxApi, DxPath, DxUtils, InstanceTypeRequest}
import dx.core.ir.RunSpec._
import dx.core.ir.{
  ExecutableKind,
  ExecutableKindNative,
  ExecutableType,
  InstanceTypeSelection,
  RuntimeRequirement,
  Value
}
import dx.core.languages.wdl.{DxRuntimeHint, IrToWdlValueBindings, Runtime, WdlUtils}

import scala.util.Try
import wdlTools.eval.WdlValues._
import wdlTools.eval.{Eval, EvalException, Meta}
import wdlTools.syntax.WdlVersion
import wdlTools.types.WdlTypes._
import wdlTools.types.{TypedAbstractSyntax => TAT}

import scala.util.matching.Regex

/**
  * A unification of WDL Runtime and Hints, with version-specific support.
  */
object RuntimeTranslator {
  // This is the parent key for the object containing all DNAnexus-specific
  // options in the hints section (2.0 and later)
  val DxKey = "dnanexus"
  // dx-specific keys that are used in meta
  val ExecutableTypeKey = "type"
  val ExecutableId = "id"
  val ExecutableTypeNative = "native"
  val AppName = "name"
  val AppletProject = "project"
  val AppletPath = "path"

  case object Access extends DxRuntimeHint(Some("dx_access"), "access", Vector(T_Object))
  case object App extends DxRuntimeHint(Some("dx_app"), "app", Vector(T_String, T_Object))
  case object IgnoreReuse
      extends DxRuntimeHint(Some("dx_ignore_reuse"), "ignore_reuse", Vector(T_Boolean))
  case object Restart extends DxRuntimeHint(Some("dx_restart"), "restart", Vector(T_Int, T_Object))
  // TODO: this is an input file hint
//  case object Stream
//      extends DxRuntimeHint(Some("dx_stream"), "stream", Vector(T_Boolean))
  case object Timeout
      extends DxRuntimeHint(Some("dx_timeout"), "timeout", Vector(T_String, T_Object))
  // TODO: case object Regions
  /**
    * This key is used in the restart object value to represent "*"
    */
  val AllKey = "All"
}

case class RuntimeTranslator(wdlVersion: WdlVersion,
                             runtimeSection: Option[TAT.RuntimeSection],
                             hintsSection: Option[TAT.MetaSection],
                             metaSection: Option[TAT.MetaSection],
                             defaultAttrs: Map[String, Value],
                             evaluator: Eval,
                             dxApi: DxApi = DxApi.get) {
  private lazy val runtime =
    Runtime(wdlVersion,
            runtimeSection,
            hintsSection,
            evaluator,
            defaultAttrs = Some(IrToWdlValueBindings(defaultAttrs)))
  private lazy val meta: Meta = Meta.create(wdlVersion, metaSection)

  def translate(id: String, wdlType: Option[T] = None): Option[Value] = {
    try {
      (runtime.get(id), wdlType) match {
        case (Some(value), None)    => Some(WdlUtils.toIRValue(value))
        case (Some(value), Some(t)) => Some(WdlUtils.toIRValue(value, t))
        case other                  => throw new Exception(s"invalid value ${other}")
      }
    } catch {
      case _: EvalException =>
        // the value is an expression that requires evaluation at runtime
        None
    }
  }

  def translateExecutableKind: Option[ExecutableKind] = {
    def kindFromId(id: String): ExecutableKind = {
      val (executableType, _) = DxUtils.parseExecutableId(id)
      ExecutableKindNative(ExecutableType.withNameIgnoreCase(executableType), Some(id))
    }
    runtime
      .getDxHint(RuntimeTranslator.App)
      .map {
        case V_String(id) =>
          try {
            kindFromId(id)
          } catch {
            case _: IllegalArgumentException =>
              if (id.startsWith("/")) {
                ExecutableKindNative(ExecutableType.Applet, path = Some(id))
              } else {
                ExecutableKindNative(ExecutableType.App, name = Some(id))
              }
          }
        case V_Object(fields) =>
          def getStringField(name: String): Option[String] = {
            fields.get(name) match {
              case Some(V_String(s)) => Some(s)
              case None              => None
              case other             => throw new Exception(s"Invalid ${name} ${other}")
            }
          }
          val id = getStringField(RuntimeTranslator.ExecutableId)
          val name = getStringField(RuntimeTranslator.AppName)
          val project = getStringField(RuntimeTranslator.AppletProject)
          val path = getStringField(RuntimeTranslator.AppletPath)
          getStringField(RuntimeTranslator.ExecutableTypeKey) match {
            case Some(executableType) =>
              ExecutableKindNative(
                  ExecutableType.withNameIgnoreCase(executableType),
                  id,
                  name,
                  project,
                  path
              )

            case None if name.isDefined =>
              ExecutableKindNative(ExecutableType.App, id, name)
            case None if project.isDefined || path.isDefined =>
              ExecutableKindNative(ExecutableType.Applet, id, project = project, path = path)
            case None if id.isDefined => kindFromId(id.get)
            case _ =>
              throw new Exception("Not enough information to determine native app(let)")
          }
        case other => throw new Exception(s"invalid app value ${other}")
      }
      .orElse {
        (meta.get(RuntimeTranslator.ExecutableTypeKey),
         meta.get(RuntimeTranslator.ExecutableId)) match {
          case (Some(V_String(RuntimeTranslator.ExecutableTypeNative)), Some(V_String(id))) =>
            Some(kindFromId(id))
          case _ => None
        }
      }
  }

  def translateInstanceType(
      resolution: InstanceTypeSelection.InstanceTypeSelection
  ): InstanceType = {
    runtime.safeParseInstanceType match {
      case Some(InstanceTypeRequest.empty) => DefaultInstanceType
      case Some(req: InstanceTypeRequest)
          if resolution == InstanceTypeSelection.Dynamic && req.dxInstanceType.isDefined =>
        StaticInstanceType(InstanceTypeRequest(dxInstanceType = req.dxInstanceType))
      case Some(_) if resolution == InstanceTypeSelection.Dynamic => DynamicInstanceType
      case Some(req: InstanceTypeRequest)                         => StaticInstanceType(req)
      case None                                                   => DynamicInstanceType
    }
  }

  def translateContainer: ContainerImage = {
    if (!runtime.containerDefined) {
      return NoImage
    }
    Try {
      // Prefer dxfile image if one is available, otherwise an
      // external Docker image if available
      val (dxUris, otherUris) = runtime.container.partition(_.startsWith(DxPath.DxUriPrefix))
      dxUris.headOption
        .map(uri => DxFileDockerImage(uri, dxApi.resolveFile(uri)))
        .orElse(otherUris.headOption.map(NetworkDockerImage))

      // If runtime.containerDefined but runtime.container is None,
      // image has to be dynamically evaluated at runtime
    }.toOption.flatten.getOrElse(DynamicDockerImage)
  }

  private def unwrapString(value: V): String = {
    value match {
      case V_String(s) => s
      case _           => throw new Exception(s"Expected V_Array(V_String), got ${value}")
    }
  }

  private def unwrapStringArray(value: V): Vector[String] = {
    value match {
      case V_Array(a) => a.map(unwrapString)
      case _          => throw new Exception(s"Expected V_String, got ${value}")
    }
  }

  private def unwrapBoolean(value: V): Boolean = {
    value match {
      case V_Boolean(b) => b
      case _            => throw new Exception(s"Expected V_Boolean, got ${value}")
    }
  }

  private def unwrapInt(value: V): Long = {
    value match {
      case V_Int(i) => i
      case _        => throw new Exception(s"Expected V_Int, got ${value}")
    }
  }

  def translateAccess: Option[AccessRequirement] = {
    runtime.getDxHint(RuntimeTranslator.Access).map {
      case V_Object(fields) =>
        AccessRequirement(
            network = fields.get("network").map(unwrapStringArray).getOrElse(Vector.empty),
            project = fields.get("project").map(unwrapString),
            allProjects = fields.get("allProjects").map(unwrapString),
            developer = fields.get("developer").map(unwrapBoolean),
            projectCreation = fields.get("projectCreation").map(unwrapBoolean)
        )
      case other => throw new Exception(s"invalid access value ${other}")
    }
  }

  def translateIgnoreReuse: Option[IgnoreReuseRequirement] = {
    runtime.getDxHint(RuntimeTranslator.IgnoreReuse).map {
      case V_Boolean(ignoreReuse) => IgnoreReuseRequirement(ignoreReuse)
      case other                  => throw new Exception(s"invalid ignoreReuse value ${other}")
    }
  }

  def translateRestart: Option[RestartRequirement] = {
    runtime.getDxHint(RuntimeTranslator.Restart).map {
      case V_Int(n) => RestartRequirement(default = Some(n))
      case V_Object(fields) =>
        RestartRequirement(
            fields.get("max").map(unwrapInt),
            fields.get("default").map(unwrapInt),
            fields
              .get("errors")
              .map {
                case V_Object(errFields) =>
                  errFields.map {
                    case (s, V_Int(i)) => s -> i
                    case other =>
                      throw new Exception(s"Invalid restart map entry ${other}")
                  }
                case _ => throw new Exception("Invalid restart map")
              }
              .getOrElse(Map.empty)
        )
      case other => throw new Exception(s"invalid restart value ${other}")
    }
  }

  private val durationRegexp = s"^(?:(\\d+)D)?(?:(\\d+)H)?(?:(\\d+)M)?".r
  // val durationFields = Vector("days", "hours", "minutes")

  private def parseDuration(duration: String): TimeoutRequirement = {
    durationRegexp.findFirstMatchIn(duration) match {
      case Some(result: Regex.Match) =>
        def group(i: Int): Option[Long] = {
          result.group(i) match {
            case null      => None
            case s: String => Some(s.toLong)
          }
        }
        TimeoutRequirement(group(1), group(2), group(3))
      case _ => throw new Exception(s"Invalid ISO Duration ${duration}")
    }
  }

  def translateTimeout: Option[TimeoutRequirement] = {
    runtime.getDxHint(RuntimeTranslator.Timeout).map {
      case V_String(s) => parseDuration(s)
      case V_Object(fields) =>
        TimeoutRequirement(
            fields.get("days").map(unwrapInt),
            fields.get("hours").map(unwrapInt),
            fields.get("minutes").map(unwrapInt)
        )
      case other => throw new Exception(s"invalid timeout value ${other}")
    }
  }

  def translateRequirements: Vector[RuntimeRequirement] = {
    Vector(
        translateAccess,
        translateIgnoreReuse,
        translateRestart,
        translateTimeout
    ).flatten
  }
}
