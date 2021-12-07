package dx.compiler

import java.time.{LocalDateTime, ZoneId}
import java.time.format.DateTimeFormatter
import dx.api.{
  DxApi,
  DxApplet,
  DxAppletDescribe,
  DxDataObject,
  DxFindDataObjects,
  DxFindDataObjectsConstraints,
  DxObjectDescribe,
  DxProject,
  DxWorkflow,
  DxWorkflowDescribe,
  Field
}
import dx.core.Constants
import dx.core.ir.Bundle
import dx.util.{JsUtils, Logger}

trait DxExecutableInfo {
  val dataObj: DxDataObject
  val desc: Option[DxObjectDescribe]
  val checksum: Option[String]
  val createdDate: Option[LocalDateTime]
  def name: String
}

object DxExecutableDirectory {
  val ArchivedTag = "archived"
}

/**
  * Takes a snapshot of the platform target path before the build starts. Provides methods for
  * querying and updating the directory.
  * @param bundle an IR bundle
  * @param project the project to search
  * @param folder the folder to search
  * @param projectWideReuse whether to allow project-wide reuse
  * @param dxApi the Dx API
  * @param logger the logger
  */
case class DxExecutableDirectory(bundle: Bundle,
                                 project: DxProject,
                                 folder: String,
                                 projectWideReuse: Boolean = false,
                                 dxApi: DxApi = DxApi.get,
                                 logger: Logger = Logger.get) {

  // a list of all dx:workflow and dx:applet names used in this WDL workflow
  private val allExecutableNames: Set[String] = bundle.allCallables.keySet
  private val dxFind = DxFindDataObjects(dxApi)

  /**
    * Information about a Dx data object.
    * @param dataObj the actual data object
    * @param dxDesc the object description
    * @param checksum the object checksum
    * @param createdDate created date
    */
  case class DxExecutableWithDesc(dataObj: DxDataObject,
                                  dxDesc: DxObjectDescribe,
                                  checksum: Option[String],
                                  createdDate: Option[LocalDateTime] = None)
      extends DxExecutableInfo {
    def name: String = dxDesc.name
    val desc: Option[DxObjectDescribe] = Some(dxDesc)
  }

  private def findExecutables(
      folder: Option[String] = None,
      recurse: Boolean = false
  ): Vector[(DxDataObject, DxObjectDescribe)] = {
    val constraints = DxFindDataObjectsConstraints(
        project = Some(project),
        folder = folder,
        recurse = recurse,
        tags = Set(Constants.CompilerTag),
        names = allExecutableNames
    )
    val dxObjectsInFolder = dxFind
      .query(
          constraints,
          withInputOutputSpec = false,
          defaultFields = true,
          extraFields = Set(Field.Tags, Field.Properties, Field.Details)
      )
      .filter {
        case (_: DxApplet, desc: DxAppletDescribe)
            if !desc.tags.exists(_.contains(DxExecutableDirectory.ArchivedTag)) =>
          true
        case (_: DxWorkflow, desc: DxWorkflowDescribe)
            if !desc.tags.exists(_.contains(DxExecutableDirectory.ArchivedTag)) =>
          true
        case _ => false
      }
      .toVector
    logger.trace(s"Found ${dxObjectsInFolder.size} executables in ${project.id}:${folder}")
    dxObjectsInFolder
  }

  private def getChecksum(desc: DxObjectDescribe): Option[String] = {
    desc.details
      .flatMap(_.asJsObject.fields.get(Constants.Checksum))
      .map(JsUtils.getString(_))
      .orElse(desc.properties.flatMap(_.get(Constants.ChecksumPropertyDeprecated)))
  }

  /*
   * Instead of looking up applets/workflows one by one, perform a bulk lookup, and find all the
   * objects in the target directory. Setup an easy to use map with information on each name.
   *
   * findDataObjects can be an expensive call, both on the server and client sides. We limit it by
   * filtering on the CHECKSUM property, which is attached only to generated applets and workflows.
   * This runs the risk of missing cases where an applet name is already in use by a regular
   * dnanexus applet/workflow.
   */
  private lazy val initialExecDir: Map[String, Vector[DxExecutableInfo]] = {
    findExecutables(Some(folder))
      .map {
        case (dxObj, desc) =>
          val creationDate = new java.util.Date(desc.created)
          val creationTime: LocalDateTime =
            LocalDateTime.ofInstant(creationDate.toInstant, ZoneId.systemDefault())
          // checksum is stored in details, but used to be stored as a property, so
          // look in both places
          val checksum = getChecksum(desc)
          DxExecutableWithDesc(dxObj, desc, checksum, Some(creationTime))
      }
      .groupBy(_.name)
  }
  private var execDir: Option[Map[String, Vector[DxExecutableInfo]]] = None

  // A map from checksum to dx:executable, across the entire project. It allows reusing executables
  // across the entire project, at the cost of a potentially expensive API call. It is not clear
  // this is useful to the majority of users, so it is gated by the [projectWideReuse] flag.
  private lazy val projectWideExecDir: Map[String, Vector[DxExecutableInfo]] = {
    if (projectWideReuse) {
      // Scan the entire project for dx:workflows and dx:applets that we already created, and may be
      // reused, instead of recompiling. This could be expensive. We limit it by filtering on the
      // CHECKSUM property, which is attached only to generated applets and workflows. The maximal
      // number of replies is (by default) 1000 so we may miss matches when we search. The cost
      // would be creating a dx:executable again, which is acceptable.
      logger.trace(s"Querying for executables in project ${project}")
      val executables = findExecutables(recurse = true)
        .flatMap {
          case (obj, desc) =>
            getChecksum(desc).map(checksum => DxExecutableWithDesc(obj, desc, Some(checksum)))
        }
        .groupBy(_.checksum.get)
      logger.trace(s"Found ${executables.size} executables")
      executables
    } else {
      Map.empty
    }
  }

  /**
    * Gets information about all executables with the given name.
    * @param name the executable name
    * @return
    */
  def lookup(name: String): Vector[DxExecutableInfo] = {
    execDir.getOrElse(initialExecDir).getOrElse(name, Vector.empty)
  }

  /**
    * Searches for an executable with a specific checksum anywhere in the project. In case of
    * checksum collision (i.e. multiple results for the same checksum), returns only the executable
    * that starts with the name we are looking for.
    * @param name the executable name to look up
    * @param digest the executable's checksum
    * @return
    */
  def lookupInProject(name: String, digest: String): Option[DxExecutableInfo] = {
    projectWideExecDir
      .get(digest)
      .flatMap(checksumMatches => checksumMatches.find(_.name.startsWith(name)))
  }

  case class DxExecutableInserted(
      name: String,
      dataObj: DxDataObject,
      checksum: Option[String],
      createdDate: Option[LocalDateTime] = Some(LocalDateTime.now)
  ) extends DxExecutableInfo {
    override val desc: Option[DxObjectDescribe] = None
  }

  /**
    * Insert an executable into the directory.
    * @param name the executable name
    * @param dxExec the data object
    * @param digest the checksum
    */
  def insert(name: String, dxExec: DxDataObject, digest: String): Unit = {
    val info = DxExecutableInserted(name, dxExec, Some(digest))
    execDir = execDir.getOrElse(initialExecDir) match {
      case d if d.contains(name) => Some(d + (name -> (d(name) :+ info)))
      case d                     => Some(d + (name -> Vector(info)))
    }
  }

  private lazy val dateFormatter = DateTimeFormatter.ofPattern("EE MMM dd kk:mm:ss yyyy")
  private var folders: Set[String] = Set.empty

  /**
    * Creates a folder if it does not already exist.
    */
  private def ensureFolder(fullPath: String): Unit = {
    if (!folders.contains(fullPath)) {
      project.newFolder(fullPath, parents = true)
      folders += fullPath
    }
  }

  /**
    * Moves an executable into an archive directory. For example if the applet is /A/B/C/GLnexus,
    * moves it to /A/B/C/.archive/GLnexus {date}/GLnexus. The archived executable is tagged with
    * "archived" so that it can be filtered out when searching for existing applets.
    * TODO: need to update the execInfo in the directory
    * @param execInfo the object to archive
    */
  def archive(execInfo: DxExecutableInfo): Unit = {
    logger.trace(s"Archiving ${execInfo.name} ${execInfo.dataObj.id}")
    // tag the object
    dxApi.addTags(execInfo.dataObj, Vector(DxExecutableDirectory.ArchivedTag))
    // move the object to the new location
    val archiveName = execInfo.createdDate match {
      case Some(dt) => s"${execInfo.name} ${dt.format(dateFormatter)}"
      case None     => execInfo.name
    }
    val destFolder = s"${folder}/.archive/${archiveName}/"
    ensureFolder(destFolder)
    project.moveObjects(Vector(execInfo.dataObj), destFolder)
  }

  /**
    * Archive multiple executables.
    * @param execInfos the executables to archive
    */
  def archive(execInfos: Vector[DxExecutableInfo]): Unit = {
    execInfos.foreach(archive)
  }

  /**
    * Remove executables from the project and update the directory.
    * TODO: need to remove the execInfos in the directory
    * @param execInfos the executables to remove
    */
  def remove(execInfos: Vector[DxExecutableInfo]): Unit = {
    val objs = execInfos.map(_.dataObj)
    logger.trace(s"Removing executables ${objs.map(_.id)}")
    project.removeObjects(objs, force = true)
  }
}
