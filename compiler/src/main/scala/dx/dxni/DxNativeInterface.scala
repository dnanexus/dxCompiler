package dx.dxni

import dx.api._
import dx.core.getVersion
import dx.core.Constants
import dx.core.languages.Language.Language
import dx.dxni.wdl.WdlDxNativeInterfaceFactory
import dx.util.{FileSourceResolver, Logger}

trait NativeInterfaceGenerator {
  def generate(apps: Vector[DxApp] = Vector.empty,
               applets: Vector[DxApplet] = Vector.empty,
               headerLines: Vector[String]): Vector[String]
}

trait NativeInterfaceGeneratorFactory {
  def create(language: Language): Option[NativeInterfaceGenerator]
}

case class DxNativeInterface(fileResolver: FileSourceResolver = FileSourceResolver.get,
                             dxApi: DxApi = DxApi.get,
                             logger: Logger = Logger.get) {

  private val generatorFactories = Vector(
      WdlDxNativeInterfaceFactory(fileResolver = fileResolver, dxApi = dxApi, logger = logger)
  )

  private def getGenerator(language: Language): NativeInterfaceGenerator = {
    generatorFactories
      .collectFirst { factory =>
        factory.create(language) match {
          case Some(support) => support
        }
      }
      .getOrElse(
          throw new Exception(s"Language ${language} is not supported")
      )
  }

  private def getApplet(dxProject: DxProject, path: String): DxApplet = {
    path match {
      case id if id.startsWith("applet-") => dxApi.applet(id)
      case _ =>
        dxApi.resolveDataObject(path, Some(dxProject)) match {
          case applet: DxApplet => applet
          case _                => throw new Exception(s"DxNI only supports apps and applets")
        }
    }
  }

  private def searchApplets(dxProject: DxProject,
                            folder: String,
                            recursive: Boolean): Vector[DxApplet] = {
    val (nativeApplets, dxcApplets): (Map[DxApplet, DxAppletDescribe],
                                      Map[DxApplet, DxAppletDescribe]) =
      DxFindDataObjects(dxApi)
        .apply(Some(dxProject),
               Some(folder),
               recursive,
               classRestriction = Some("applet"),
               withInputOutputSpec = true,
               extraFields = Set(Field.Tags))
        .collect {
          case (applet: DxApplet, desc: DxAppletDescribe) => (applet, desc)
        }
        .partition {
          // ignore any applets with the compiler tag set - it indicates an applet
          // that was compiled with dxCompiler (and thus not "native")
          case (_: DxApplet, desc: DxAppletDescribe)
              if !desc.tags.get.contains(Constants.CompilerTag) =>
            true
          case (_: DxApplet, _: DxAppletDescribe) => false
        }

    if (dxcApplets.nonEmpty) {
      val dxcAppletNames = dxcApplets.map {
        case (_, desc: DxAppletDescribe) => desc.name
      }
      logger.warning(
          s"Applets ${dxcAppletNames.mkString(", ")} were ignored by dxni because they were compiled with dxCompiler"
      )
    }
    if (nativeApplets.isEmpty) {
      logger.warning(s"Found no applets in ${dxProject.id}/${folder}")
    } else {
      logger.trace(s"Found ${nativeApplets.size} applets in ${dxProject.id}/${folder}")
    }
    nativeApplets.map {
      case (applet: DxApplet, _) => applet
    }.toVector
  }

  private def searchApps: Vector[DxApp] = {
    val apps: Vector[DxApp] = DxFindApps(dxApi)
      .apply(published = Some(true), withInputOutputSpec = true)
    if (apps.isEmpty) {
      logger.warning(s"Found no DX global apps")
    } else {
      logger.trace(s"Found ${apps.size} DX global apps")
    }
    apps
  }

  private val appsHeader = Vector(
      s"This file was generated by the Dx Native Interface (DxNI) tool ${getVersion}.",
      "These are interfaces to apps."
  )

  /**
    * Generate only apps.
    */
  def apply(language: Language, appId: Option[String]): (Vector[DxApp], Vector[String]) = {
    val generator = getGenerator(language)
    val apps: Vector[DxApp] = appId.map(id => Vector(DxApp(id)(dxApi))).getOrElse(searchApps)
    if (apps.nonEmpty) {
      (apps, generator.generate(apps, headerLines = appsHeader))
    } else {
      (Vector.empty, Vector.empty)
    }
  }

  private def appletsHeader(dxProject: DxProject, path: String) = Vector(
      s"This file was generated by the Dx Native Interface (DxNI) tool ${getVersion}.",
      s"project name = ${dxProject.describe().name}",
      s"project ID = ${dxProject.id}",
      s"path = ${path}"
  )

  def apply(language: Language,
            dxProject: DxProject,
            folder: Option[String] = None,
            path: Option[String] = None,
            applet: Option[DxApplet] = None,
            recursive: Boolean = false,
            includeApps: Boolean = true): (Vector[DxApp], Vector[DxApplet], Vector[String]) = {
    val generator = getGenerator(language)
    val apps = if (includeApps) {
      searchApps
    } else {
      Vector.empty
    }
    val (applets: Vector[DxApplet], search) = (folder, path, applet) match {
      case (Some(folder), None, None) => (searchApplets(dxProject, folder, recursive), folder)
      case (None, Some(path), None)   => (Vector(getApplet(dxProject, path)), path)
      case (None, None, Some(applet)) => (Vector(applet), applet.id)
      case _                          => throw new Exception("must specify exactly one of (folder, path)")
    }
    if (apps.nonEmpty || applets.nonEmpty) {
      (apps, applets, generator.generate(apps, applets, appletsHeader(dxProject, search)))
    } else {
      (Vector.empty, Vector.empty, Vector.empty)
    }
  }
}
