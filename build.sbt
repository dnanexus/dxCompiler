import Merging.customMergeStrategy
import sbt.Keys._
import sbtassembly.AssemblyPlugin.autoImport.{assemblyMergeStrategy, _}
import com.typesafe.config._

name := "dxc"

ThisBuild / organization := "com.dnanexus"
ThisBuild / scalaVersion := "2.13.2"
ThisBuild / developers := List(
    Developer("jdidion", "jdidion", "jdidion@dnanexus.com", url("https://github.com/dnanexus")),
    Developer("commandlinegirl",
              "commandlinegirl",
              "azalcman@dnanexus.com",
              url("https://github.com/dnanexus")),
    Developer("mhrvol", "mhrvol", "mhrvol-cf@dnanexus.com", url("https://github.com/dnanexus"))
)
ThisBuild / homepage := Some(url("https://github.com/dnanexus/dxCompiler"))
ThisBuild / scmInfo := Some(
    ScmInfo(url("https://github.com/dnanexus/dxCompiler"), "git@github.com:dnanexus/dxCompiler.git")
)
ThisBuild / licenses += ("Apache-2.0", url("http://www.apache.org/licenses/LICENSE-2.0"))

// PROJECTS

lazy val root = project.in(file("."))
lazy val global = root
  .settings(
      settings,
      skip in publish := true
  )
  .disablePlugins(AssemblyPlugin)
  .aggregate(
      core,
      compiler,
      executorCommon,
      executorWdl,
      executorCwl
  )

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
          dependencies.cwlScala
      ),
      assemblyJarName in assembly := "dxCompiler.jar",
      assemblyOutputPath in assembly := file("applet_resources/dxCompiler.jar")
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
      assemblyJarName in assembly := "dxExecutorWdl.jar",
      assemblyOutputPath in assembly := file("applet_resources/WDL/resources/dxExecutorWdl.jar")
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
      assemblyJarName in assembly := "dxExecutorCwl.jar",
      assemblyOutputPath in assembly := file("applet_resources/CWL/resources/dxExecutorCwl.jar")
  )
  .dependsOn(core, executorCommon)

// DEPENDENCIES

val githubDxScalaResolver = Resolver.githubPackages("dnanexus", "dxScala")
val githubCwlScalaResolver = Resolver.githubPackages("dnanexus", "cwlScala")
val githubWdlToolsResolver = Resolver.githubPackages("dnanexus-rnd", "wdlTools")
val githubDxCompilerResolver = Resolver.githubPackages("dnanexus", "dxCompiler")

lazy val dependencies =
  new {
    val dxCommonVersion = "0.6.1-SNAPSHOT"
    val dxApiVersion = "0.7.1-SNAPSHOT"
    val dxFileAccessProtocolsVersion = "0.4.1"
    val wdlToolsVersion = "0.14.4-SNAPSHOT"
    val cwlScalaVersion = "0.5.1-SNAPSHOT"
    val typesafeVersion = "1.3.3"
    val sprayVersion = "1.3.5"
    val scalatestVersion = "3.1.1"
    val logbackVersion = "1.2.3"

    val dxCommon = "com.dnanexus" % "dxcommon" % dxCommonVersion
    val dxApi = "com.dnanexus" % "dxapi" % dxApiVersion
    val dxFileAccessProtocols = "com.dnanexus" % "dxfileaccessprotocols" % dxFileAccessProtocolsVersion
    val wdlTools = "com.dnanexus" % "wdltools" % wdlToolsVersion
    val cwlScala = "com.dnanexus" % "cwlscala" % cwlScalaVersion
    val typesafe = "com.typesafe" % "config" % typesafeVersion
    val spray = "io.spray" %% "spray-json" % sprayVersion
    val logback = "ch.qos.logback" % "logback-classic" % logbackVersion
    val scalatest = "org.scalatest" % "scalatest_2.13" % scalatestVersion
  }

lazy val commonDependencies = Seq(
    dependencies.dxCommon,
    dependencies.dxApi,
    dependencies.dxFileAccessProtocols,
    dependencies.logback,
    dependencies.spray,
    dependencies.scalatest % Test
)

// SETTINGS

lazy val settings = Seq(
    scalacOptions ++= compilerOptions,
    // javac
    javacOptions ++= Seq("-Xlint:deprecation"),
    // reduce the maximum number of errors shown by the Scala compiler
    maxErrors := 20,
    // scalafmt
    scalafmtConfig := root.base / ".scalafmt.conf",
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
    Test / testOptions += Tests.Argument("-oF"),
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
    "-Xlint:nullary-override",
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
    logLevel in assembly := Level.Info,
    // comment out this line to enable tests in assembly
    test in assembly := {},
    assemblyMergeStrategy in assembly := customMergeStrategy.value
)
