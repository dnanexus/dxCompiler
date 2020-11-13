package dx.api

import spray.json._
import dx.util.Enum

// Extra fields for describe
object Field extends Enum {
  type Field = Value
  val Access, Analysis, App, Applet, ArchivalState, AvailableInstanceTypes, BillTo, Categories,
      Created, Description, Details, DeveloperNotes, Executable, ExecutableName, Folder, Id,
      IgnoreReuse, Input, Inputs, InputSpec, InstanceType, Modified, Name, Output, Outputs,
      OutputSpec, ParentJob, Parts, PricingModelsByRegion, Project, Properties, Region, RunSpec,
      Size, Stages, Summary, Tags, Title, Types, Version = Value
}

trait DxObjectDescribe {
  val id: String
  val name: String
  val created: Long
  val modified: Long
  val properties: Option[Map[String, String]]
  val details: Option[JsValue]

  def getCreationDate: java.util.Date = new java.util.Date(created)
}

trait DxObject {
  def id: String

  def describe(fields: Set[Field.Value] = Set.empty): DxObjectDescribe
}

object DxObject {
  def parseJsonProperties(props: JsValue): Map[String, String] = {
    props.asJsObject.fields.map {
      case (k, JsString(v)) => k -> v
      case (_, _) =>
        throw new Exception(s"malformed JSON properties ${props}")
    }
  }

  def maybeSpecifyProject(project: Option[DxProject]): Map[String, JsValue] = {
    project match {
      case None =>
        // we don't know the project.
        Map.empty
      case Some(p) =>
        // We know the project, this makes the search more efficient.
        Map("project" -> JsString(p.id))
    }
  }

  def requestFields(fields: Set[Field.Value]): JsValue = {
    val fieldStrings = fields.map {
      case Field.Access                 => "access"
      case Field.Analysis               => "analysis"
      case Field.App                    => "app"
      case Field.Applet                 => "applet"
      case Field.ArchivalState          => "archivalState"
      case Field.AvailableInstanceTypes => "availableInstanceTypes"
      case Field.BillTo                 => "billTo"
      case Field.Categories             => "categories"
      case Field.Created                => "created"
      case Field.Description            => "description"
      case Field.DeveloperNotes         => "developerNotes"
      case Field.Details                => "details"
      case Field.Executable             => "executable"
      case Field.ExecutableName         => "executableName"
      case Field.Folder                 => "folder"
      case Field.Id                     => "id"
      case Field.IgnoreReuse            => "ignoreReuse"
      case Field.Input                  => "input"
      case Field.Inputs                 => "inputs"
      case Field.InputSpec              => "inputSpec"
      case Field.InstanceType           => "instanceType"
      case Field.Modified               => "modified"
      case Field.Name                   => "name"
      case Field.Output                 => "output"
      case Field.Outputs                => "outputs"
      case Field.OutputSpec             => "outputSpec"
      case Field.ParentJob              => "parentJob"
      case Field.Parts                  => "parts"
      case Field.PricingModelsByRegion  => "pricingModelsByRegion"
      case Field.Project                => "project"
      case Field.Properties             => "properties"
      case Field.Region                 => "region"
      case Field.RunSpec                => "runSpec"
      case Field.Size                   => "size"
      case Field.Stages                 => "stages"
      case Field.Summary                => "summary"
      case Field.Tags                   => "tags"
      case Field.Title                  => "title"
      case Field.Types                  => "types"
      case Field.Version                => "version"
    }.toVector
    val m: Map[String, JsValue] = fieldStrings.map { x =>
      x -> JsTrue
    }.toMap
    JsObject(m)
  }
}

trait DxDataObject extends DxObject

// Objects that can be run on the platform. These are apps, applets, and workflows.
trait DxExecutable extends DxObject

// Actual executions on the platform. There are jobs and analyses
trait DxExecution extends DxObject {
  val id: String
  val project: Option[DxProject]
}

// DxDataObject that caches its description
abstract class CachingDxObject[T <: DxObjectDescribe] extends DxObject {
  private var cachedDesc: Option[T] = None

  def hasCachedDesc: Boolean = cachedDesc.nonEmpty

  def cacheDescribe(desc: DxObjectDescribe): Unit = {
    cachedDesc = Some(desc.asInstanceOf[T])
  }

  def describe(fields: Set[Field.Value] = Set.empty): T = {
    if (cachedDesc.isEmpty) {
      cachedDesc = Some(describeNoCache(fields))
    } else {
      // TODO: check that all `fields` are present in desc
    }
    cachedDesc.get
  }

  def describeNoCache(fields: Set[Field.Value]): T
}
