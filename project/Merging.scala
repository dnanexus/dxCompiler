// Copied verbatim from the Cromwell project
// https://github.com/broadinstitute/cromwell/blob/develop/project/Merging.scala
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtassembly.{MergeStrategy, PathList}

case class OnErrorMergeStrategy(strategy: MergeStrategy, fallback: MergeStrategy)
    extends MergeStrategy {
  override def name: String = s"OnError(${strategy.name},${fallback.name})"

  def apply(tempDir: File, path: String, files: Seq[File]): Either[String, Seq[(File, String)]] = {
    // TODO: log failure message
    strategy.apply(tempDir, path, files) match {
      case Right(success) => Right(success)
      case Left(_)        => fallback.apply(tempDir, path, files)
    }
  }
}

object Merging {
  val customMergeStrategy: Def.Initialize[String => MergeStrategy] = Def.setting {
    case PathList(ps @ _*) if ps.last == "project.properties" =>
      // Merge/Filter project.properties files from Google jars that otherwise collide at merge time.
      MergeStrategy.filterDistinctLines
    case PathList(ps @ _*) if ps.last == "logback.xml" =>
      MergeStrategy.first
    // AWS SDK v2 configuration files - can be discarded
    case PathList(ps @ _*)
        if Set("codegen.config",
               "service-2.json",
               "waiters-2.json",
               "customization.config",
               "examples-1.json",
               "paginators-1.json").contains(ps.last) =>
      MergeStrategy.discard
    case x @ PathList("META-INF", path @ _*) =>
      path map {
        _.toLowerCase
      } match {
        case "spring.tooling" :: _ =>
          MergeStrategy.discard
        case "io.netty.versions.properties" :: Nil =>
          MergeStrategy.first
        case "maven" :: "com.google.guava" :: _ =>
          MergeStrategy.first
        case _ :+ ("asm-license.txt" | "module-info.class" | "overview.html" |
            "cobertura.properties") =>
          MergeStrategy.discard
        case _ :+ ("generated-layer.xml" | "namedservices.index" | "java.lang.runnable") =>
          MergeStrategy.discard
        case _ =>
          val oldStrategy = (assembly / assemblyMergeStrategy).value
          oldStrategy(x)
      }
    case x @ PathList("OSGI-INF", path @ _*) =>
      path map {
        _.toLowerCase
      } match {
        case "l10n" :: "bundle.properties" :: Nil =>
          MergeStrategy.concat
        case _ =>
          val oldStrategy = (assembly / assemblyMergeStrategy).value
          oldStrategy(x)
      }
    case PathList(
        "asm-license.txt" | "module-info.class" | "overview.html" | "cobertura.properties"
        ) =>
      MergeStrategy.discard
    case PathList("mime.types") =>
      MergeStrategy.last
    case x =>
      val oldStrategy = (assembly / assemblyMergeStrategy).value
      OnErrorMergeStrategy(oldStrategy(x), MergeStrategy.first)
  }
}
