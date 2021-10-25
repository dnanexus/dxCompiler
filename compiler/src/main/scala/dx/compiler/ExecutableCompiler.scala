package dx.compiler

import dx.api.{DxApi, DxConstraint, DxConstraintBool, DxConstraintString, DxIOSpec}
import dx.core.Constants
import dx.core.ir.Value._
import dx.core.ir.{
  Callable,
  DxName,
  Parameter,
  ParameterAttribute,
  ParameterLinkSerializer,
  Type,
  TypeSerde,
  Value,
  ValueSerde
}
import dx.translator.CallableAttributes.{
  DescriptionAttribute,
  DetailsAttribute,
  PropertiesAttribute,
  SummaryAttribute,
  TagsAttribute,
  TitleAttribute,
  TypesAttribute
}
import dx.translator.{Extras, ParameterAttributes}
import dx.util.FileSourceResolver
import spray.json._

object ExecutableCompiler {
  // these parameters are used for applets and workflows that are generated
  // with useManifests=true.
  val InputManifestParameter: Parameter =
    Parameter(Constants.InputManifest, Type.TOptional(Type.THash))
  val InputManfestFilesParameter: Parameter =
    Parameter(Constants.InputManifestFiles, Type.TArray(Type.TFile))
  val InputLinksParameter: Parameter =
    Parameter(Constants.InputLinks, Type.TOptional(Type.THash))
  val WorkflowInputManifestParameter: Parameter =
    Parameter(Constants.WorkflowInputManifest, Type.TOptional(Type.THash))
  val WorkflowInputManfestFilesParameter: Parameter =
    Parameter(Constants.WorkflowInputManifestFiles, Type.TArray(Type.TFile))
  val WorkflowInputLinksParameter: Parameter =
    Parameter(Constants.WorkflowInputLinks, Type.TOptional(Type.THash))
  val OutputIdParameter: Parameter =
    Parameter(Constants.OutputId, Type.TString)
  val CallNameParameter: Parameter =
    Parameter(Constants.CallName, Type.TOptional(Type.TString))
  val OutputManifestParameter: Parameter =
    Parameter(Constants.OutputManifest, Type.TFile)
  val OverridesParameter: Parameter =
    Parameter(Constants.Overrides, Type.TOptional(Type.THash))
}

class ExecutableCompiler(extras: Option[Extras],
                         parameterLinkSerializer: ParameterLinkSerializer,
                         complexPathValues: Boolean,
                         fileResolver: FileSourceResolver,
                         dxApi: DxApi = DxApi.get) {

  private def typeConstraintToNative(constraint: DxConstraint): JsValue = {
    constraint match {
      case DxConstraintString(s) => JsString(s)
      case DxConstraintBool(oper, array) =>
        JsObject(Map(oper.name -> JsArray(array.constraints.map(typeConstraintToNative))))
      case other =>
        throw new Exception(s"invalid type constraint ${other}")
    }
  }

  // Create the IO Attributes
  private def parameterAttributesToNative(attrs: Vector[ParameterAttribute],
                                          dxType: Type,
                                          excludeAttributes: Set[String]): Map[String, JsValue] = {
    attrs
      .flatMap {
        case ParameterAttributes.GroupAttribute(text) =>
          Some(DxIOSpec.Group -> JsString(text))
        case ParameterAttributes.HelpAttribute(text) =>
          Some(DxIOSpec.Help -> JsString(text))
        case ParameterAttributes.LabelAttribute(text) =>
          Some(DxIOSpec.Label -> JsString(text))
        case ParameterAttributes.PatternsAttribute(patterns) =>
          patterns match {
            case ParameterAttributes.PatternsArray(array) =>
              Some(DxIOSpec.Patterns -> JsArray(array.map(JsString(_))))
            // If we have the alternative patterns object, extrac the values, if any at all
            case ParameterAttributes.PatternsObject(name, klass, tags) =>
              Vector(
                  if (name.isEmpty) None else Some("name" -> JsArray(name.map(JsString(_)))),
                  if (tags.isEmpty) None else Some("tag" -> JsArray(tags.map(JsString(_)))),
                  klass.map("class" -> JsString(_))
              ).flatten match {
                case Vector() => None
                case v        => Some(DxIOSpec.Patterns -> JsObject(v.toMap))
              }
          }
        case ParameterAttributes.ChoicesAttribute(choices) =>
          Some(DxIOSpec.Choices -> JsArray(choices.collect {
            case ParameterAttributes.SimpleChoice(VString(value))  => JsString(value)
            case ParameterAttributes.SimpleChoice(VInt(value))     => JsNumber(value)
            case ParameterAttributes.SimpleChoice(VFloat(value))   => JsNumber(value)
            case ParameterAttributes.SimpleChoice(VBoolean(value)) => JsBoolean(value)
            case ParameterAttributes.FileChoice(value, name) => {
              // TODO: support project and record choices
              val dxLink = dxApi.resolveFile(value).asJson
              if (name.isDefined) {
                JsObject(Map("name" -> JsString(name.get), "value" -> dxLink))
              } else {
                dxLink
              }
            }
            // TODO: ParameterAttributes.DirectoryChoice
          }))
        case ParameterAttributes.SuggestionsAttribute(suggestions) =>
          Some(DxIOSpec.Suggestions -> JsArray(suggestions.collect {
            case ParameterAttributes.SimpleSuggestion(VString(value))  => JsString(value)
            case ParameterAttributes.SimpleSuggestion(VInt(value))     => JsNumber(value)
            case ParameterAttributes.SimpleSuggestion(VFloat(value))   => JsNumber(value)
            case ParameterAttributes.SimpleSuggestion(VBoolean(value)) => JsBoolean(value)
            case ParameterAttributes.FileSuggestion(value, name, project, path) => {
              // TODO: support project and record suggestions
              val dxLink: Option[JsValue] = value match {
                case Some(str) => Some(dxApi.resolveFile(str).asJson)
                case None      => None
              }
              if (name.isDefined || project.isDefined || path.isDefined) {
                val attrs: Map[String, JsValue] = Vector(
                    if (dxLink.isDefined) Some("value" -> dxLink.get) else None,
                    if (name.isDefined) Some("name" -> JsString(name.get)) else None,
                    if (project.isDefined) Some("project" -> JsString(project.get)) else None,
                    if (path.isDefined) Some("path" -> JsString(path.get)) else None
                ).flatten.toMap
                JsObject(attrs)
              } else if (dxLink.isDefined) {
                dxLink.get
              } else {
                throw new Exception(
                    "Either 'value' or 'project' + 'path' must be defined for suggestions"
                )
              }
            }
            // TODO: ParameterAttributes.DirectorySuggestion
          }))
        case ParameterAttributes.TypeAttribute(constraint) =>
          Some(DxIOSpec.Type -> typeConstraintToNative(constraint))
        case ParameterAttributes.DefaultAttribute(value) =>
          // The default was specified in parameter_meta and was not specified in the
          // parameter specification
          Some(
              DxIOSpec.Default -> ValueSerde.serializeWithType(value,
                                                               dxType,
                                                               fileResolver = Some(fileResolver))
          )
        case _ => None
      }
      .toMap
      .filterNot {
        case (key, _) => excludeAttributes.contains(key)
      }
  }

  private def optionalToNative(optional: Boolean): Map[String, JsValue] = {
    if (optional) {
      Map("optional" -> JsTrue)
    } else {
      Map.empty[String, JsValue]
    }
  }

  private val InputOnlyKeys: Set[String] = {
    Set(DxIOSpec.Default, DxIOSpec.Choices, DxIOSpec.Suggestions)
  }

  /**
    * Converts an IR Paramter to a native input/output spec.
    * - For primitive types, and arrays of such types, we can map directly
    *   to the equivalent dx types. For example,
    *     Int  -> int
    *     Array[String] -> array:string
    * - Arrays can be empty, which is why they are always marked "optional".
    *   This notifies the platform runtime system not to throw an exception
    *   for an empty input/output array.
    * - Ragged arrays, maps, and objects, cannot be mapped in such a trivial way.
    *   These are called "Complex Types", or "Complex". They are handled
    *   by passing a JSON structure and a vector of dx:files.
    * - For parameters with defaults, we always set them to optional regardless
    *   of their WdlType.
    * @param parameter the input Parameter
    * @return the DNAnexus inputDesc
    */
  protected def inputParameterToNative(parameter: Parameter): Vector[JsObject] = {
    val defaultValues: Map[DxName, JsValue] = parameter.defaultValue match {
      case Some(value) =>
        parameterLinkSerializer.createFields(parameter.name, parameter.dxType, value).toMap
      case None => Map.empty
    }

    def defaultValueToNative(name: DxName): Map[String, JsValue] = {
      defaultValues.get(name) match {
        case Some(jsv) => Map(DxIOSpec.Default -> jsv)
        case None      => Map.empty
      }
    }

    val excludeAttributeNames: Set[String] = if (defaultValues.contains(parameter.name)) {
      Set(DxIOSpec.Default)
    } else {
      Set.empty
    }
    val attributes = defaultValueToNative(parameter.name) ++
      parameterAttributesToNative(parameter.attributes, parameter.dxType, excludeAttributeNames)
    val (nativeType, optional) = TypeSerde.toNative(parameter.dxType, !complexPathValues)
    // TODO: I don't think the parameter should always be set to optional if it has a default
    val paramSpec = JsObject(
        Map(DxIOSpec.Name -> JsString(parameter.name.encoded),
            DxIOSpec.Class -> JsString(nativeType))
          ++ attributes
          ++ optionalToNative(optional || attributes.contains(DxIOSpec.Default))
    )
    if (nativeType == "hash") {
      // A JSON structure passed as a hash, and a vector of platform files
      val filesName = parameter.name.addSuffix(Constants.FlatFilesSuffix)
      Vector(
          paramSpec,
          JsObject(
              Map(
                  DxIOSpec.Name -> JsString(filesName.encoded),
                  DxIOSpec.Class -> JsString("array:file"),
                  DxIOSpec.Optional -> JsTrue
              )
                ++ defaultValueToNative(filesName)
              // some attributes don't makes sense for the files input
                ++ attributes.filterNot(x => InputOnlyKeys.contains(x._1))
          )
      )
    } else {
      Vector(paramSpec)
    }
  }

  /**
    * Similar to inputParameterToNative, but output parameters don't allow some
    * fields: default, suggestions, choices.
    * @param parameter the output Parameter
    * @return the DNAnexus outputDesc
    */
  protected def outputParameterToNative(parameter: Parameter): Vector[JsObject] = {
    val attributes =
      parameterAttributesToNative(parameter.attributes, parameter.dxType, InputOnlyKeys)
    val (nativeType, optional) = TypeSerde.toNative(parameter.dxType, !complexPathValues)
    val paramSpec = JsObject(
        Map(DxIOSpec.Name -> JsString(parameter.name.encoded),
            DxIOSpec.Class -> JsString(nativeType))
          ++ optionalToNative(optional || parameter.defaultValue.isDefined)
          ++ attributes
    )
    if (nativeType == "hash") {
      // A JSON structure passed as a hash, and a vector of platform files
      val filesName = parameter.name.addSuffix(Constants.FlatFilesSuffix)
      Vector(
          paramSpec,
          JsObject(
              Map(
                  DxIOSpec.Name -> JsString(filesName.encoded),
                  DxIOSpec.Class -> JsString("array:file"),
                  DxIOSpec.Optional -> JsTrue
              )
                ++ attributes
          )
      )
    } else {
      Vector(paramSpec)
    }
  }

  // Match everything up to the first period; truncate after 50 characters.
  private val MaxSummaryLength = 50
  private lazy val firstLineRegex = s"^([^.]{1,${MaxSummaryLength}}).*".r

  private def summaryToNative(summary: Option[String],
                              description: Option[String]): Map[String, JsValue] = {
    (summary, description) match {
      case (Some(text), _) if text.nonEmpty                 => Map("summary" -> JsString(text))
      case (_, Some(firstLineRegex(line))) if line.nonEmpty =>
        // Default 'summary' to be the first line of 'description'
        val descSummary = if (line.length() == MaxSummaryLength && !line.endsWith(".")) {
          line + "..."
        } else {
          line
        }
        Map("summary" -> JsString(descSummary))
      case _ => Map.empty
    }
  }

  // Merge initial description & extended description
  private def descriptionToNative(description: Option[String],
                                  extendedDescription: Option[String]): Map[String, JsValue] = {
    (description, extendedDescription) match {
      case (Some(d), Some(e)) => Map("description" -> JsString(s"${d}\n\n${e}"))
      case (Some(d), _)       => Map("description" -> JsString(d))
      case (_, Some(e))       => Map("description" -> JsString(e))
      case (_, _)             => Map.empty
    }
  }

  private def whatsNewToNative(whatsNew: Option[JsValue]): Map[String, JsValue] = {
    whatsNew match {
      case Some(JsArray(array)) =>
        // If whatsNew is in array format, convert it to a string
        val changelog = array
          .map {
            case JsObject(fields) =>
              val formattedFields = fields
                .map {
                  case ("version", JsString(value)) => Some("version" -> value)
                  case ("changes", JsArray(array)) =>
                    Some("changes",
                         array
                           .map {
                             case JsString(item) => s"* ${item}"
                             case other =>
                               throw new Exception(s"Invalid change list item: ${other}")
                           }
                           .mkString("\n"))
                  case _ => None
                }
                .flatten
                .toMap
              s"### Version ${formattedFields("version")}\n${formattedFields("changes")}"
            case other => throw new Exception(s"Invalid whatsNew ${other}")
          }
          .mkString("\n")
        Map("whatsNew" -> JsString(s"## Changelog\n${changelog}"))
      case _ => Map.empty
    }
  }

  // Convert the applet meta to JSON, and overlay details from task-specific extras
  protected def callableAttributesToNative(
      callable: Callable,
      defaultTags: Set[String],
      extendedDescription: Option[String]
  ): (Map[String, JsValue], Map[String, JsValue]) = {
    val metaDefaults = Map(
        "title" -> JsString(callable.name),
        "tags" -> JsArray((defaultTags ++ callable.tags).map(JsString(_)).toVector),
        "properties" -> JsObject(callable.properties.map {
          case (name, value) => name -> JsString(value)
        })
        // These are currently ignored because they only apply to apps
        //"version" -> JsString("0.0.1"),
        //"openSource" -> JsBoolean(false),
    )
    val meta = callable.attributes.collect {
      case TitleAttribute(text)  => "title" -> JsString(text)
      case TypesAttribute(array) => "types" -> JsArray(array.map(JsString(_)))
      case TagsAttribute(array)  =>
        // merge default and user-specified tags
        "tags" -> JsArray((array.toSet ++ defaultTags ++ callable.tags).map(JsString(_)).toVector)
      case PropertiesAttribute(props) =>
        "properties" -> JsObject((props ++ callable.properties).map {
          case (k, v) => k -> JsString(v)
        })
    }.toMap
    // String attributes that need special handling
    val meta2: Map[String, String] = callable.attributes.collect {
      case SummaryAttribute(text)     => "summary" -> text
      case DescriptionAttribute(text) => "description" -> text
    }.toMap
    // default summary to the first line of the description
    val summary =
      summaryToNative(meta2.get("summary"), meta2.get("description"))
    // Merge initial description & extended description
    val description = descriptionToNative(meta2.get("description"), extendedDescription)
    // extract the details to return separately
    val metaDetails = callable.attributes
      .collectFirst {
        case DetailsAttribute(details) =>
          details.map {
            case (k, v) => k -> ValueSerde.serialize(v)
          }
      }
      .getOrElse(Map.empty)
    // get the whatsNew section from the details
    val whatsNew = whatsNewToNative(metaDetails.get("whatsNew"))
    (metaDefaults ++ meta ++ summary ++ description, metaDetails ++ whatsNew)
  }

  protected def delayWorkspaceDestructionToNative: Map[String, JsValue] = {
    if (extras.flatMap(_.delayWorkspaceDestruction).getOrElse(false)) {
      Map(Constants.DelayWorkspaceDestruction -> JsTrue)
    } else {
      Map.empty
    }
  }

  protected def fileDependenciesFromParam(param: Parameter): Set[String] = {
    val default = param.defaultValue match {
      case Some(f: Value.VFile) => Vector(f.uri)
      case _                    => Vector.empty
    }
    val attrs = param.attributes
      .flatMap(attr =>
        attr match {
          case ParameterAttributes.DefaultAttribute(f: Value.VFile) => Vector(f.uri)
          case ParameterAttributes.ChoicesAttribute(choices) =>
            choices
              .collect {
                case ParameterAttributes.SimpleChoice(f: Value.VFile) => f.uri
                case ParameterAttributes.FileChoice(str, _)           => str
              }
          case ParameterAttributes.SuggestionsAttribute(suggestions) =>
            suggestions
              .collect {
                case ParameterAttributes.SimpleSuggestion(f: Value.VFile)   => f.uri
                case ParameterAttributes.FileSuggestion(Some(str), _, _, _) => str
              }
          case _ => Vector.empty
        }
      )

    (default ++ attrs).filter(s => s.contains("://")).toSet
  }

  // TODO: Use templates for Markdown dependency report

  /**
    * Markdown helper method for 1st-level list item
    * @param s text string
    * @return Markdown list item of text string
    */
  protected def listMd(s: String): String = {
    s"\n* ${s}"
  }

  /**
    * Markdown helper method for 2nd-level indented list item, as code
    * @param s text string
    * @return Markdown 2nd-level list item of text string, as code
    */
  protected def listMd2(s: String): String = {
    s"\n    * `${s}`"
  }
}
