import Merging.customMergeStrategy
import sbt.Keys._
import com.typesafe.config._

name := "dxc"

ThisBuild / organization := "com.dnanexus"
ThisBuild / scalaVersion := "2.13.7"
ThisBuild / developers := List(
    Developer("commandlinegirl",
              "Ola Zalcman",
              "azalcman@dnanexus.com",
              url("https://github.com/dnanexus")),
    Developer("Gvaihir", "Gvaihir", "aogrodnikov@dnanexus.com", url("https://github.com/dnanexus")),
    Developer("mhrvol", "Marek Hrvol", "mhrvol@dnanexus.com", url("https://github.com/dnanexus")),
    Developer("r-i-v-a",
              "Riva Nathans",
              "rnathans@dnanexus.com",
              url("https://github.com/dnanexus")),
    Developer("YuxinShi0423", "Yuxin Shi", "yshi@dnanexus.com", url("https://github.com/dnanexus")),
)
ThisBuild / homepage := Some(url("https://github.com/dnanexus/dxCompiler"))
ThisBuild / scmInfo := Some(
    ScmInfo(url("https://github.com/dnanexus/dxCompiler"), "git@github.com:dnanexus/dxCompiler.git")
)
ThisBuild / licenses += ("Apache-2.0", url("http://www.apache.org/licenses/LICENSE-2.0"))

// PROJECTS

lazy val root = project
  .in(file("."))
  .settings(
      settings,
      publish / skip := true
  )
  .aggregate(
      core,
      compiler,
      executorCommon,
      executorWdl,
      executorCwl
  )
  .disablePlugins(AssemblyPlugin)

val dxCompilerVersion: String = {
  val confPath = s"core/src/main/resources/application.conf"
  val conf = ConfigFactory.parseFile(new File(confPath)).resolve()
  conf.getString(s"dxCompilerCore.version")
}

val core = project
  .in(file("core"))
  .settings(
      name := "dxCompilerCore",
      version := dxCompilerVersion,
      settings,
      libraryDependencies ++= commonDependencies ++ Seq(
          dependencies.typesafe,
          dependencies.wdlTools,
          dependencies.cwlScala
      )
  )
  .disablePlugins(AssemblyPlugin)

val compiler = project
  .in(file("compiler"))
  .settings(
      name := "dxCompiler",
      version := dxCompilerVersion,
      settings,
      assemblySettings,
      libraryDependencies ++= commonDependencies ++ Seq(
          dependencies.typesafe,
          dependencies.wdlTools,
          dependencies.cwlScala,
          dependencies.dxYaml
      ),
      assembly / assemblyJarName := "dxCompiler.jar",
      assembly / assemblyOutputPath := file("applet_resources/dxCompiler.jar")
  )
  .dependsOn(core)

val executorCommon = project
  .in(file("executorCommon"))
  .settings(
      name := "dxExecutorCommon",
      version := dxCompilerVersion,
      settings,
      libraryDependencies ++= commonDependencies ++ Seq()
  )
  .disablePlugins(AssemblyPlugin)
  .dependsOn(core)

val executorWdl = project
  .in(file("executorWdl"))
  .settings(
      name := "dxExecutorWdl",
      version := dxCompilerVersion,
      settings,
      assemblySettings,
      libraryDependencies ++= commonDependencies ++ Seq(
          dependencies.typesafe,
          dependencies.wdlTools
      ),
      assembly / assemblyJarName := "dxExecutorWdl.jar",
      assembly / assemblyOutputPath := file("applet_resources/WDL/resources/dxExecutorWdl.jar")
  )
  .dependsOn(core, executorCommon)

val executorCwl = project
  .in(file("executorCwl"))
  .settings(
      name := "dxExecutorCwl",
      version := dxCompilerVersion,
      settings,
      assemblySettings,
      libraryDependencies ++= commonDependencies ++ Seq(
          dependencies.typesafe,
          dependencies.cwlScala
      ),
      assembly / assemblyJarName := "dxExecutorCwl.jar",
      assembly / assemblyOutputPath := file("applet_resources/CWL/resources/dxExecutorCwl.jar")
  )
  .dependsOn(core, executorCommon)

// DEPENDENCIES

val githubDxScalaResolver = Resolver.githubPackages("dnanexus", "dxScala")
val githubCwlScalaResolver = Resolver.githubPackages("dnanexus", "cwlScala")
val githubWdlToolsResolver = Resolver.githubPackages("dnanexus", "wdlTools")
val githubDxCompilerResolver = Resolver.githubPackages("dnanexus", "dxCompiler")

lazy val dependencies =
  new {
    val dxCommonVersion = "0.11.5"
    val dxApiVersion = "0.13.8"
    val dxFileAccessProtocolsVersion = "0.5.6"
    val dxYamlVersion = "0.1.1"
    val cwlScalaVersion = "0.8.4"
    val wdlToolsVersion = "0.17.16"
    val typesafeVersion = "1.4.1"
    val sprayVersion = "1.3.6"
    val scalatestVersion = "3.2.9"
    val logbackVersion = "1.2.10"
    val mockitoVersion = "3.2.10.0"

    val dxCommon = "com.dnanexus" % "dxcommon" % dxCommonVersion
    val dxApi = "com.dnanexus" % "dxapi" % dxApiVersion
    val dxFileAccessProtocols = "com.dnanexus" % "dxfileaccessprotocols" % dxFileAccessProtocolsVersion
    val dxYaml = "com.dnanexus" % "dxyaml" % dxYamlVersion
    val wdlTools = "com.dnanexus" % "wdltools" % wdlToolsVersion
    val cwlScala = "com.dnanexus" % "cwlscala" % cwlScalaVersion
    val typesafe = "com.typesafe" % "config" % typesafeVersion
    val spray = "io.spray" %% "spray-json" % sprayVersion
    val logback = "ch.qos.logback" % "logback-classic" % logbackVersion
    val scalatest = "org.scalatest" % "scalatest_2.13" % scalatestVersion
    val mockito = "org.scalatestplus" %% "mockito-3-4" % mockitoVersion % "test"
  }

lazy val commonDependencies = Seq(
    dependencies.dxCommon,
    dependencies.dxApi,
    dependencies.dxFileAccessProtocols,
    dependencies.logback,
    dependencies.spray,
    dependencies.scalatest % Test,
    dependencies.mockito
)

// SETTINGS

// exclude tests tagged 'linuxOnly' unless we're on a Linux OS
val isLinux = System.getProperty("os.name").toLowerCase.contains("linux")
val scalatestArgs = if (isLinux) {
  Seq("-oFK")
} else {
  Seq("-oFK", "-l", "linuxOnly")
}
lazy val settings = Seq(
    scalacOptions ++= compilerOptions,
    // javac
    javacOptions ++= Seq("-Xlint:deprecation"),
    // reduce the maximum number of errors shown by the Scala compiler
    maxErrors := 20,
    // scalafmt
    scalafmtConfig := file(".") / ".scalafmt.conf",
    // Publishing
    // disable publish with scala version, otherwise artifact name will include scala version
    // e.g dxScala_2.11
    crossPaths := false,
    // snapshot artifact resolvers
    resolvers ++= Seq(githubDxScalaResolver,
                      githubWdlToolsResolver,
                      githubCwlScalaResolver,
                      githubDxCompilerResolver),
    // add sonatype repository settings
    // snapshot versions publish to GitHub packages repository
    // release versions publish to sonatype staging repository
    publishTo := Some(
        if (isSnapshot.value) {
          githubDxCompilerResolver
        } else {
          Opts.resolver.sonatypeStaging
        }
    ),
    githubOwner := "dnanexus",
    githubRepository := "dxCompiler",
    publishMavenStyle := true,
    // Tests
    // If an exception is thrown during tests, show the full
    // stack trace, by adding the "-oF" option to the list.
    // Ignore cancelled tests (-oK)
    Test / testOptions += Tests.Argument(scalatestArgs: _*),
    Test / parallelExecution := false
    // Coverage
    //
    // sbt clean coverage test
    // sbt coverageReport
    //coverageEnabled := true
    // To turn it off do:
    // sbt coverageOff
    // Ignore code parts that cannot be checked in the unit
    // test environment
    //coverageExcludedPackages := "dxScala.Main"
)

// Show deprecation warnings
val compilerOptions = Seq(
    "-unchecked",
    "-deprecation",
    "-feature",
    "-explaintypes",
    "-encoding",
    "UTF-8",
    "-Xlint:constant",
    "-Xlint:delayedinit-select",
    "-Xlint:doc-detached",
    "-Xlint:inaccessible",
    "-Xlint:infer-any",
    "-Xlint:nullary-unit",
    "-Xlint:option-implicit",
    "-Xlint:package-object-classes",
    "-Xlint:poly-implicit-overload",
    "-Xlint:private-shadow",
    "-Xlint:stars-align",
    "-Xlint:type-parameter-shadow",
    "-Ywarn-dead-code",
    "-Ywarn-unused:implicits",
    "-Ywarn-unused:privates",
    "-Ywarn-unused:locals",
    "-Ywarn-unused:imports", // warns about every unused import on every command.
    "-Xfatal-warnings" // makes those warnings fatal.
)

// Assembly
lazy val assemblySettings = Seq(
    // comment out this line to enable tests in assembly
    assembly / test := {},
    assembly / assemblyMergeStrategy := customMergeStrategy.value
)
