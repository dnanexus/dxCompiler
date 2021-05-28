package dxCompiler

import java.nio.file.{Files, Path, Paths}
import com.typesafe.config.{Config, ConfigFactory}
import dx.api._
import dx.compiler.{Compiler, ExecutableTree}
import dx.core.getVersion
import dx.core.CliUtils._
import dx.core.io.{DxWorkerPaths, StreamFiles}
import dx.core.ir.Bundle
import dx.core.languages.Language
import dx.core.Constants
import dx.core.languages.wdl.WdlOptions
import dx.dxni.DxNativeInterface
import dx.translator.{Extras, TranslatorFactory}
import dx.util.protocols.DxFileAccessProtocol
import dx.util.{Enum, FileSourceResolver, FileUtils, Logger, TraceLevel}
import spray.json.{JsNull, JsValue}
import wdlTools.types.TypeCheckingRegime

/**
  * Compiler CLI.
  */
object Main {
  private val DefaultRuntimeTraceLevel: Int = TraceLevel.Verbose

  /**
    * OptionSpec that parses the string argument(s) to -language.
    */
  private case class LanguageOptionSpec() extends OptionSpec {

    /**
      * Parses a language argument. Accepts the following:
      *  'cwl'
      *  'cwlv1.2'
      *  'cwl 1.2'
      *  'draft2' -> WDL draft2
      */
    override def parseValues(name: String, values: Vector[String], curValue: Option[Opt]): Opt = {
      if (curValue.nonEmpty) {
        throw OptionParseException(s"Option ${name} specified multiple times")
      }
      val language = values match {
        case Vector(language) =>
          Language.parse(language)
        case Vector(language, version) =>
          Language.parse(version, Some(language))
        case _ =>
          throw OptionParseException(s"Unexpected value ${values} to option ${name}")
      }
      SingleValueOption[Language.Language](language)
    }
  }

  private val CommonOptions: InternalOptions = Map(
      "destination" -> StringOptionSpec.one,
      "force" -> FlagOptionSpec.default,
      "f" -> FlagOptionSpec.default.copy(alias = Some("force")),
      "overwrite" -> FlagOptionSpec.default.copy(alias = Some("force")),
      "folder" -> StringOptionSpec.one,
      "project" -> StringOptionSpec.one,
      "language" -> LanguageOptionSpec()
  )

  /**
    * Initialize objects used commonly by multiple commands.
    * - creates a FileSourceResolver that looks for local files in any configured -imports
    *   directories and has a DxFileAccessProtocol
    * - initializes a Logger
    * @param options parsed options
    * @return (FileSourceResolver, Logger)
    */
  private def initCommon(options: Options): (FileSourceResolver, Logger) = {
    val logger = initLogger(options)
    val imports: Vector[Path] = options.getList[Path]("imports")
    val fileResolver = FileSourceResolver.create(
        imports,
        Vector(DxFileAccessProtocol()),
        logger
    )
    FileSourceResolver.set(fileResolver)
    (fileResolver, logger)
  }

  // compile

  object CompilerAction extends Enum {
    type CompilerAction = Value
    val Compile, Config, DxNI, Version, Describe = Value
  }

  object CompilerMode extends Enum {
    type CompilerMode = Value
    val IR, NativeWithoutRuntimeAsset, All = Value
  }

  private case class CompilerModeOptionSpec()
      extends SingleValueOptionSpec[CompilerMode.CompilerMode](
          choices = CompilerMode.values.toVector
      ) {
    override def parseValue(value: String): CompilerMode.CompilerMode =
      CompilerMode.withNameIgnoreCase(value)
  }

  // Tree printer types for the execTree option
  object ExecTreeFormat extends Enum {
    type ExecTreeFormat = Value
    val Json, Pretty = Value
  }

  private case class ExecTreeFormatOptionSpec()
      extends SingleValueOptionSpec[ExecTreeFormat.ExecTreeFormat](
          choices = ExecTreeFormat.values.toVector
      ) {
    override def parseValue(value: String): ExecTreeFormat.ExecTreeFormat =
      ExecTreeFormat.withNameIgnoreCase(value)
  }

  private object StreamFilesOptionSpec
      extends SingleValueOptionSpec[StreamFiles.StreamFiles](choices = StreamFiles.values.toVector) {
    override def parseValue(value: String): StreamFiles.StreamFiles = {
      StreamFiles.withNameIgnoreCase(value)
    }
  }

//  private object WdlRegimeOptionSpec
//      extends SingleValueOptionSpec[TypeCheckingRegime.TypeCheckingRegime](
//          choices = TypeCheckingRegime.values.toVector
//      ) {
//    override def parseValue(value: String): TypeCheckingRegime.TypeCheckingRegime = {
//      TypeCheckingRegime.withNameIgnoreCase(value)
//    }
//  }

  private def CompileOptions: InternalOptions = Map(
      "archive" -> FlagOptionSpec.default,
      "compileMode" -> CompilerModeOptionSpec(),
      "defaults" -> PathOptionSpec.mustExist,
      "execTree" -> ExecTreeFormatOptionSpec(),
      "extras" -> PathOptionSpec.mustExist,
      "inputs" -> PathOptionSpec.listMustExist,
      "input" -> PathOptionSpec.listMustExist.copy(alias = Some("inputs")),
      "locked" -> FlagOptionSpec.default,
      "leaveWorkflowsOpen" -> FlagOptionSpec.default,
      "imports" -> PathOptionSpec.listMustExist,
      "p" -> PathOptionSpec.listMustExist.copy(alias = Some("imports")),
      "projectWideReuse" -> FlagOptionSpec.default,
      "reorg" -> FlagOptionSpec.default,
      "runtimeDebugLevel" -> IntOptionSpec.one.copy(choices = Vector(0, 1, 2)),
      "separateOutputs" -> FlagOptionSpec.default,
      "streamFiles" -> StreamFilesOptionSpec,
      "streamAllFiles" -> FlagOptionSpec.default,
      "scatterChunkSize" -> IntOptionSpec.one,
      "useManifests" -> FlagOptionSpec.default,
      "waitOnUpload" -> FlagOptionSpec.default
      //"wdlMode" -> WdlRegimeOptionSpec
  )

  private val DeprecatedCompileOptions = Set(
      "fatalValidationWarnings",
      "input",
      "streamAllFiles"
  )

  private def resolvePath(
      dxApi: DxApi,
      project: Option[String],
      folder: Option[String],
      path: Option[String] = None,
      create: Boolean = false
  ): (DxProject, Either[String, DxDataObject]) = {
    val projectId = project.getOrElse({
      val projectId = dxApi.currentProjectId.get
      Logger.get.warning(s"Project is unspecified...using currently selected project ${projectId}")
      projectId
    })
    val dxProject =
      try {
        dxApi.resolveProject(projectId)
      } catch {
        case t: Throwable =>
          throw new Exception(
              s"Could not find project ${project}, you probably need to be logged into the platform",
              t
          )
      }
    val folderOrPath = (folder, path) match {
      case (Some(f), None) =>
        // Validate the folder.
        // TODO: check for folder existance rather than listing the contents, which could
        //   be very large.
        if (create) {
          dxProject.newFolder(f, parents = true)
        } else {
          dxProject.listFolder(f)
        }
        Left(f)
      case (None, Some(p)) =>
        // validate the file
        val dataObj = dxApi.resolveDataObject(p, Some(dxProject))
        Right(dataObj)
      case (None, None) =>
        Left("/")
      case _ =>
        throw OptionParseException("must specify only one of (folder, path)")
    }
    Logger.get.trace(s"""|project ID: ${dxProject.id}
                         |path: ${folderOrPath}""".stripMargin)
    (dxProject, folderOrPath)
  }

  private def resolveOrCreatePath(dxApi: DxApi,
                                  project: String,
                                  folder: String): (DxProject, String) = {
    resolvePath(dxApi, Some(project), Some(folder), create = true) match {
      case (dxProject, Left(folder)) => (dxProject, folder)
      case _                         => throw new Exception("expected folder")
    }
  }

  // There are three possible syntaxes:
  //    project-id:/folder
  //    project-id:
  //    /folder
  private def getDestination(dxApi: DxApi, options: Options): (DxProject, String) = {
    val destinationOpt: Option[String] = options.getValue[String]("destination")
    val folderOpt: Option[String] = options.getValue[String]("folder")
    val projectOpt: Option[String] = options.getValue[String]("project")
    val destRegexp = "(.+):(.*)".r
    val (project, folder) = (destinationOpt, projectOpt, folderOpt) match {
      case (Some(destRegexp(project, folder)), _, _) if folder.startsWith("/") =>
        (project, folder)
      case (Some(destRegexp(project, emptyFolder)), _, Some(folder)) if emptyFolder.trim.isEmpty =>
        (project, folder)
      case (Some(destRegexp(project, emptyFolder)), _, None) if emptyFolder.trim.isEmpty =>
        (project, "/")
      case (Some(destRegexp(_, folder)), _, None) =>
        throw OptionParseException(s"Invalid folder <${folder}>")
      case (Some(folder), Some(project), _) if folder.startsWith("/") =>
        (project, folder)
      case (Some(folder), None, _) if folder.startsWith("/") && dxApi.currentProjectId.isDefined =>
        val project = dxApi.currentProjectId.get
        Logger.get.warning(s"Project is unspecified...using currently selected project ${project}")
        (project, folder)
      case (Some(other), _, _) =>
        throw OptionParseException(s"Invalid destination <${other}>")
      case (None, Some(project), Some(folder)) =>
        (project, folder)
      case (None, Some(project), None) =>
        (project, "/")
      case (None, None, Some(folder)) =>
        val project = dxApi.currentProjectId.get
        Logger.get.warning(s"Project is unspecified...using currently selected project ${project}")
        (project, folder)
      case (None, None, None) =>
        val project = dxApi.currentProjectId.get
        Logger.get.warning(s"Project is unspecified...using currently selected project ${project}")
        (project, "/")
      case _ =>
        throw OptionParseException("Project is unspecified")
    }
    resolveOrCreatePath(dxApi, project, folder)
  }

  sealed trait CompilerSuccessfulTermination extends SuccessfulTermination {
    def compilerAction: CompilerAction.CompilerAction
  }

  sealed trait SuccessfulCompile extends CompilerSuccessfulTermination {
    override val compilerAction: CompilerAction.CompilerAction = CompilerAction.Compile

    def compilerMode: CompilerMode.CompilerMode
  }

  case class SuccessfulCompileIR(bundle: Bundle) extends SuccessfulCompile {
    override val compilerMode: CompilerMode.CompilerMode = CompilerMode.IR

    override val message: String = "Intermediate representation"
  }

  sealed trait SuccessfulPrettyTree extends CompilerSuccessfulTermination {
    def prettyTree: String
  }

  sealed trait SuccessfulJsonTree extends CompilerSuccessfulTermination {
    def jsonTree: JsValue
  }

  sealed trait SuccessfulCompileNative extends SuccessfulCompile {
    def executableIds: Vector[String]

    override def message: String = {
      executableIds.mkString(",")
    }
  }

  case class SuccessfulCompileNativeNoTree(override val compilerMode: CompilerMode.CompilerMode,
                                           override val executableIds: Vector[String])
      extends SuccessfulCompileNative

  case class SuccessfulCompileNativeWithPrettyTree(
      override val compilerMode: CompilerMode.CompilerMode,
      override val executableIds: Vector[String],
      override val prettyTree: String
  ) extends SuccessfulCompileNative
      with SuccessfulPrettyTree

  case class SuccessfulCompileNativeWithJsonTree(
      override val compilerMode: CompilerMode.CompilerMode,
      override val executableIds: Vector[String],
      override val jsonTree: JsValue
  ) extends SuccessfulCompileNative
      with SuccessfulJsonTree

  def compile(args: Vector[String]): Termination = {
    val sourceFile: Path = args.headOption
      .map(Paths.get(_))
      .getOrElse(
          throw OptionParseException(
              "Missing required positional argument <WDL file>"
          )
      )
    val options: Options =
      try {
        parseCommandLine(args.tail, CommonOptions ++ CompileOptions, DeprecatedCompileOptions)
      } catch {
        case e: OptionParseException =>
          return BadUsageTermination("Error parsing command line options", Some(e))
      }

    val (baseFileResolver, logger) = initCommon(options)
    val dxApi = DxApi()(logger)

    val extras: Option[Extras] =
      options.getValue[Path]("extras").map(extrasPath => Extras.parse(extrasPath))
    if (extras.exists(_.customReorgAttributes.isDefined)) {
      val conflictingOpts = Set("reorg", "locked").filter(options.contains)
      if (conflictingOpts.nonEmpty) {
        throw OptionParseException(
            s"ERROR: cannot provide --reorg option when ${conflictingOpts.mkString(",")} is specified in extras."
        )
      }
    }

    val defaultScatterChunkSize: Int = options.getValue[Int]("scatterChunkSize") match {
      case None => Constants.JobPerScatterDefault
      case Some(size) =>
        if (size < 1) {
          Constants.JobPerScatterDefault
        } else if (size > Constants.JobsPerScatterLimit) {
          logger.warning(
              s"The number of jobs per scatter must be between 1-${Constants.JobsPerScatterLimit}"
          )
          Constants.JobsPerScatterLimit
        } else {
          size
        }
    }

    val compileMode: CompilerMode.CompilerMode =
      options.getValueOrElse[CompilerMode.CompilerMode]("compileMode", CompilerMode.All)

    val useManifests: Boolean = options.getFlag("useManifests")
    val locked = options.getFlag("locked") match {
      case false if useManifests =>
        logger.warning(
            "Usage of manifests is only compatible with locked workflows; setting '-locked' flag"
        )
        true
      case b => b
    }

//    val wdlOptions = options
//      .getValue[TypeCheckingRegime.TypeCheckingRegime]("wdlMode")
//      .map(regime => WdlOptions(regime))
//      .getOrElse(WdlOptions.default)
    val wdlOptions = WdlOptions.default

    if (wdlOptions.regime == TypeCheckingRegime.Lenient) {
      logger.warning(
          """You have enabled 'lenient' WDL mode, which allows the compilation of
            |tasks and workflows that contain syntax that is not compliant with the
            |WDL specification, but which may be used by some 'industry-standard'
            |workflows. This may result in execution errors. Please report any
            |issues that prevent you from using 'moderate' or 'strict' mode to the
            |author of the non-compliant WDL.""".stripMargin
      )
    }

    val translator =
      try {
        val language = options.getValue[Language.Language]("language")
        val reorg = options.getFlag("reorg")
        TranslatorFactory.createTranslator(
            sourceFile,
            language,
            wdlOptions,
            extras,
            defaultScatterChunkSize,
            locked,
            if (reorg) Some(true) else None,
            useManifests,
            baseFileResolver
        )
      } catch {
        case e: Throwable =>
          return Failure(s"Error creating translator for ${sourceFile}", exception = Some(e))
      }

    // generate IR
    val rawBundle =
      try {
        translator.apply
      } catch {
        case e: Throwable =>
          return Failure(s"Error translating ${sourceFile} to IR", exception = Some(e))
      }

    // if there are inputs they need to be translated to dx inputs
    val inputs: Vector[Path] = options.getList[Path]("inputs")
    // if there are defaults, they need to be "embedded" in the bundle
    val defaults: Option[Path] = options.getValue[Path]("defaults")
    val hasInputs = inputs.nonEmpty || defaults.nonEmpty

    // quit here if the target is IR and there are no inputs to translate
    if (!hasInputs && compileMode == CompilerMode.IR) {
      return SuccessfulCompileIR(rawBundle)
    }

    // for everything past this point, the user needs to be logged in
    if (!dxApi.isLoggedIn) {
      return Failure(s"You must be logged in to compile using mode ${compileMode}")
    }

    // a destination is only required if we are doing input translation and/or
    // compiling native apps
    val (project, folder) =
      try {
        getDestination(dxApi, options)
      } catch {
        case optEx: OptionParseException =>
          return BadUsageTermination(exception = Some(optEx))
        case ex: Throwable =>
          return Failure("Could not resolve destination", Some(ex))
      }

    val (bundle, fileResolver) = if (hasInputs) {
      val (bundleWithDefaults, fileResolver) =
        try {
          translator.translateInputs(rawBundle, inputs, defaults, project)
        } catch {
          case ex: Throwable =>
            return Failure("Error translating inputs", Some(ex))
        }
      if (compileMode == CompilerMode.IR) {
        // if we're only performing translation to IR, we can quit early
        return SuccessfulCompileIR(bundleWithDefaults)
      }
      (bundleWithDefaults, fileResolver)
    } else {
      (rawBundle, baseFileResolver)
    }

    try {
      val dxPathConfig = DxWorkerPaths.default
      val runtimeTraceLevel: Int =
        options.getValueOrElse[Int]("runtimeDebugLevel", DefaultRuntimeTraceLevel)
      val includeAsset = compileMode == CompilerMode.All
      val Vector(
          archive,
          force,
          leaveWorkflowsOpen,
          locked,
          projectWideReuse,
          separateOutputs,
          streamAllFiles,
          waitOnUpload
      ) = Vector(
          "archive",
          "force",
          "leaveWorkflowsOpen",
          "locked",
          "projectWideReuse",
          "separateOutputs",
          "streamAllFiles",
          "waitOnUpload"
      ).map(options.getFlag(_))
      val streamFiles = options.getValue[StreamFiles.StreamFiles]("streamFiles") match {
        case Some(value)            => value
        case None if streamAllFiles => StreamFiles.All
        case None                   => StreamFiles.PerFile
      }
      val compiler = Compiler(
          extras,
          dxPathConfig,
          runtimeTraceLevel,
          includeAsset,
          translator.runtimeAssetName,
          translator.runtimeJar,
          archive,
          force,
          leaveWorkflowsOpen,
          locked,
          projectWideReuse,
          separateOutputs,
          streamFiles,
          waitOnUpload,
          useManifests,
          translator.complexPathValues,
          fileResolver
      )
      val results = compiler.apply(bundle, project, folder)
      // generate the execution tree if requested
      (results.primary, options.getValue[ExecTreeFormat.ExecTreeFormat]("execTree")) match {
        case (Some(primary), Some(format)) =>
          val treeJs = primary.execTree match {
            case Some(execTree) => execTree
            case None           => ExecutableTree(results.executables).apply(primary)
          }
          format match {
            case ExecTreeFormat.Json =>
              SuccessfulCompileNativeWithJsonTree(compileMode, results.executableIds, treeJs)
            case ExecTreeFormat.Pretty =>
              SuccessfulCompileNativeWithPrettyTree(compileMode,
                                                    results.executableIds,
                                                    ExecutableTree.prettyPrint(treeJs.asJsObject))
          }
        case _ =>
          SuccessfulCompileNativeNoTree(compileMode, results.executableIds)
      }
    } catch {
      case e: Throwable =>
        Failure(exception = Some(e))
    }

  }

  // DxNI

  object AppsOption extends Enum {
    type AppsOption = Value
    val Include, Exclude, Only = Value
  }

  private def DxNIOptions: InternalOptions = Map(
      "appsOnly" -> FlagOptionSpec.default,
      "apps" -> StringOptionSpec(choices = AppsOption.names.map(_.toLowerCase).toVector),
      "path" -> StringOptionSpec.one,
      "outputFile" -> PathOptionSpec.default,
      "output" -> PathOptionSpec.default.copy(alias = Some("outputFile")),
      "o" -> PathOptionSpec.default.copy(alias = Some("outputFile")),
      "recursive" -> FlagOptionSpec.default,
      "r" -> FlagOptionSpec.default.copy(alias = Some("recursive"))
  )

  case class SuccessfulDxNI(appsOption: AppsOption.AppsOption,
                            apps: Vector[DxApp],
                            applets: Vector[DxApplet])
      extends CompilerSuccessfulTermination {
    override val compilerAction: CompilerAction.CompilerAction = CompilerAction.DxNI

    override def message: String = {
      s"Apps: ${apps.size}, Applets: ${applets.size}"
    }
  }

  def dxni(args: Vector[String]): Termination = {
    val options =
      try {
        parseCommandLine(args, CommonOptions ++ DxNIOptions)
      } catch {
        case e: OptionParseException =>
          return BadUsageTermination("Error parsing command line options", Some(e))
      }
    val (fileResolver, logger) = initCommon(options)
    val dxApi = DxApi()(logger)

    // make sure the user is logged in
    if (!dxApi.isLoggedIn) {
      return Failure(s"You must be logged in to generate stubs for native app(let)s")
    }

    val language = options.getValue[Language.Language]("language").getOrElse(Language.WdlDefault)
    val outputPath: Option[Path] = options.getValue[Path]("outputFile")
    val projectOpt = options.getValue[String]("project")
    val folderOpt = options.getValue[String]("folder")
    val pathOpt = options.getValue[String]("path")
    val pathIsAppId = pathOpt.exists { path =>
      try {
        val (objClass, _) = DxUtils.parseObjectId(path)
        objClass == "app"
      } catch {
        case _: Throwable => false
      }
    }
    // flags
    val Vector(
        appsOnly,
        force,
        recursive
    ) = Vector(
        "appsOnly",
        "force",
        "recursive"
    ).map(options.getFlag(_))
    val appsOption = if (pathIsAppId) {
      AppsOption.Only
    } else {
      options
        .getValue[String]("apps")
        .map(AppsOption.withNameIgnoreCase)
        .getOrElse(
            if (appsOnly) {
              AppsOption.Only
            } else if (Vector(projectOpt, folderOpt, pathOpt).exists(_.isDefined)) {
              AppsOption.Exclude
            } else {
              AppsOption.Include
            }
        )
    }

    def writeOutput(doc: Vector[String]): Unit = {
      val path = outputPath.map { path =>
        if (Files.exists(path)) {
          if (!force) {
            throw new Exception(
                s"Output file ${outputPath.toString} already exists, use -force to overwrite it"
            )
          }
          path.toFile.delete
        }
        path
      }
      if (path.isDefined) {
        FileUtils.writeFileContent(path.get, doc.mkString("\n"))
      } else {
        System.out.println(doc.mkString("\n"))
      }
    }

    val dxni = DxNativeInterface(fileResolver)
    if (appsOption == AppsOption.Only) {
      try {
        val (apps, lines) = dxni.apply(language, pathOpt)
        writeOutput(lines)
        SuccessfulDxNI(appsOption, apps, Vector.empty)
      } catch {
        case e: Throwable => Failure(exception = Some(e))
      }
    } else {
      val (dxProject, folderOrFile) = resolvePath(dxApi, projectOpt, folderOpt, pathOpt)
      val includeApps = appsOption match {
        case AppsOption.Include => true
        case AppsOption.Exclude => false
        case _ =>
          throw new RuntimeException(s"unexpected value for --apps ${appsOption}")
      }
      try {
        val (apps, applets, lines) = folderOrFile match {
          case Left(folder) =>
            dxni.apply(language,
                       dxProject,
                       folder = Some(folder),
                       recursive = recursive,
                       includeApps = includeApps)
          case Right(applet: DxApplet) =>
            dxni.apply(language, dxProject, applet = Some(applet), includeApps = includeApps)
          case _ =>
            throw OptionParseException(
                s"Invalid folder/path ${folderOrFile}"
            )
        }
        writeOutput(lines)
        SuccessfulDxNI(appsOption, apps, applets)
      } catch {
        case e: Throwable => Failure(exception = Some(e))
      }
    }
  }

  // describe

  sealed trait SuccessfulDescribe extends CompilerSuccessfulTermination {
    override def compilerAction: CompilerAction.CompilerAction = CompilerAction.Describe
  }

  case class SuccessfulDescribePrettyTree(override val prettyTree: String)
      extends SuccessfulDescribe
      with SuccessfulPrettyTree {
    override def message: String = prettyTree
  }

  case class SuccessfulDescribeJsonTree(override val jsonTree: JsValue)
      extends SuccessfulDescribe
      with SuccessfulJsonTree {
    override def message: String = jsonTree match {
      case JsNull => ""
      case _      => jsonTree.prettyPrint
    }
  }

  private def DescribeOptions: InternalOptions = Map(
      "pretty" -> FlagOptionSpec.default
  )

  def describe(args: Vector[String]): Termination = {
    val workflowId = args.headOption.getOrElse(
        throw OptionParseException(
            "Missing required positional argument <WDL file>"
        )
    )
    val options =
      try {
        parseCommandLine(args.tail, DescribeOptions, DeprecatedCompileOptions)
      } catch {
        case e: OptionParseException =>
          return BadUsageTermination("Error parsing command line options", Some(e))
      }
    val logger = initLogger(options)
    val dxApi = DxApi()(logger)
    // make sure the user is logged in
    if (!dxApi.isLoggedIn) {
      return Failure(s"You must be logged in to generate stubs to describe a workflow")
    }
    try {
      val wf = dxApi.workflow(workflowId)
      val execTreeJS = ExecutableTree.fromDxWorkflow(wf)
      if (options.getFlag("pretty")) {
        val prettyTree = ExecutableTree.prettyPrint(execTreeJS.asJsObject)
        SuccessfulDescribePrettyTree(prettyTree)
      } else {
        SuccessfulDescribeJsonTree(execTreeJS)
      }
    } catch {
      case e: Throwable =>
        BadUsageTermination(exception = Some(e))
    }
  }

  case class SuccessfulConfig(config: Config) extends CompilerSuccessfulTermination {
    override def compilerAction: CompilerAction.CompilerAction = CompilerAction.Config

    override def message: String = config.toString
  }

  case class SuccessfulVersion(version: String = getVersion) extends CompilerSuccessfulTermination {
    override def compilerAction: CompilerAction.CompilerAction = CompilerAction.Version

    override def message: String = version
  }

  private def dispatchCommand(args: Vector[String]): Termination = {
    if (args.isEmpty) {
      return BadUsageTermination()
    }
    val action =
      try {
        CompilerAction.withNameIgnoreCase(args.head.replaceAll("_", ""))
      } catch {
        case _: NoSuchElementException =>
          return BadUsageTermination()
      }
    try {
      action match {
        case CompilerAction.Compile  => compile(args.tail)
        case CompilerAction.Describe => describe(args.tail)
        case CompilerAction.DxNI     => dxni(args.tail)
        case CompilerAction.Config   => SuccessfulConfig(ConfigFactory.load())
        case CompilerAction.Version  => SuccessfulVersion()
      }
    } catch {
      case e: Throwable =>
        BadUsageTermination(exception = Some(e))
    }
  }

  /*
   Add the following if/when wdlMode is enabled

         -wdlMode [lenient,moderate,strict]
                             Strictness to use when parsing WDL documents. The default
                             ('moderate') will suffice for all workflows that are
                             compliant with the WDL specification, while 'lenient' will
                             enable compilation of third-party workflows with
                             non-compliant syntax.
   */
  private val usageMessage =
    s"""|java -jar dxCompiler.jar <action> <parameters> [options]
        |
        |Actions:
        |  version
        |    Prints the dxCompiler version.
        |  
        |  config
        |    Prints the current dxCompiler configuration.
        |  
        |  describe <DxWorkflow ID>
        |    Generate the JSON execution tree for a given DNAnexus workflow ID.
        |    The workflow needs to be have been previously compiled by dxCompiler.
        |    options
        |      -pretty                Print exec tree in "pretty" text format instead of JSON.
        |
        |  compile <WDL file>
        |    Compile a WDL file to a DNAnexus workflow or applet.
        |    options
        |      -archive               Archive older versions of applets.
        |      -compileMode <string>  Compilation mode - a debugging flag for internal use.
        |      -defaults <string>     JSON file with standard-formatted default values.
        |      -destination <string>    Full platform path (project:/folder).
        |      -execTree [json,pretty]    
        |                             Print a JSON representation of the workflow.
        |      -extras <string>       JSON file with extra options (see documentation).
        |      -inputs <string>       JSON file with standard-formatted input values. May be
        |                             specified multiple times. A DNAnexus JSON input file is
        |                             generated for each standard input file.
        |      -locked                Create a locked workflow. When running a locked workflow,
        |                             input values may only be specified for the top-level workflow.
        |      -leaveWorkflowsOpen    Leave created workflows open (otherwise they are closed).
        |      -p | -imports <string> Directory to search for imported WDL files. May be specified
        |                             multiple times.
        |      -projectWideReuse      Look for existing applets/workflows in the entire project
        |                             before generating new ones. The default search scope is the
        |                             target folder only.
        |      -reorg                 Reorganize workflow output files.
        |      -runtimeDebugLevel [0,1,2] 
        |                             How much debug information to write to the job log at runtime.
        |                             Log the minimum (0), intermediate (1, the default), or all 
        |                             debug information (2, for internal debugging).
        |      -separateOutputs       Store the output files of each call in a separate folder. The
        |                             default behavior is to put all outputs in the same folder.                             
        |      -streamFiles [all,none,perfile]
        |                             Whether to mount all files with dxfuse (do not use the 
        |                             download agent), to mount no files with dxfuse (only use 
        |                             download agent), or to respect the per-file settings in WDL
        |                             parameter_meta sections (default).
        |      -useManifests          Use manifest files for all workflow and applet inputs and 
        |                             outputs. Implies -locked.
        |      -waitOnUpload          Whether to wait for each file upload to complete.
        |
        |  dxni
        |    DNAnexus Native call Interface. Creates stubs for calling DNAnexus executables 
        |    (apps/applets/workflows), and stores them as WDL tasks in a local file. Enables 
        |    calling existing platform executables without modification.
        |    options:
        |      -apps [include,exclude,only]
        |                             Whether to 'include' apps, 'exclude' apps (the default), or 
        |                             'only' generate app stubs.
        |      -f | force             Delete any existing output file.
        |      -o <path>              Destination file for WDL task definitions (defaults to 
        |                             stdout).
        |      -path <string>         Name of a specific app or a path to a specific applet.
        |      -r | recursive         Search recursively for applets in the target folder.
        |
        |Common options
        |    -folder <string>         Platform folder (defaults to '/').
        |    -project <string>        Platform project (defaults to currently selected project).
        |    -language <string> [ver] Which language to use? May be WDL or CWL. You can optionally 
        |                             specify a version. Currently, WDL draft-2, 1.0, and 1.1 are
        |                             fully supported and WDL development and CWL 1.2 are partially
        |                             supported. The default is to auto-detect the language from the
        |                             source file.
        |    -quiet                   Do not print warnings or informational outputs.
        |    -verbose                 Print detailed logging.
        |    -verboseKey <module>     Print verbose output only for a specific module. May be 
        |                             specified multiple times.
        |    -logFile <path>          File to use for logging output; defaults to stderr.
        |""".stripMargin

  def main(args: Vector[String]): Unit = {
    terminate(dispatchCommand(args), usageMessage)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
