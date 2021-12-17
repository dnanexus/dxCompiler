addSbtPlugin("com.eed3si9n" % "sbt-assembly" % "1.1.0")
addSbtPlugin("com.github.sbt" % "sbt-release" % "1.1.0")
//addSbtPlugin("org.scoverage" % "sbt-scoverage" % "1.5.1")
addSbtPlugin("org.scalameta" % "sbt-scalafmt" % "2.4.3")
// sbt-github-packages plugin used to publish snapshot versions to github packages
addSbtPlugin("com.codecommit" % "sbt-github-packages" % "0.5.2")
addSbtPlugin("net.virtual-void" % "sbt-dependency-graph" % "0.10.0-RC1")

// workaround for missing static SLF4J binder for logback
libraryDependencies += "ch.qos.logback" % "logback-classic" % "1.2.8"
// only load git plugin if we're actually in a git repo
libraryDependencies ++= {
  if ((baseDirectory.value / "../.git").isDirectory)
    Seq(
        Defaults.sbtPluginExtra("com.typesafe.sbt" % "sbt-git" % "1.0.0",
                                (update / sbtBinaryVersion).value,
                                (update / scalaBinaryVersion).value)
    )
  else {
    println("sbt-git plugin not loaded")
    Seq.empty
  }
}
