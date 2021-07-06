package dx.translator

import dx.api.DxApi
import dx.core.CliUtils.UnsuccessfulTermination
import dx.util.{FileSourceResolver, JsUtils, LocalFileSource, Logger}
import dxCompiler.Main
import dxCompiler.Main.SuccessfulCompileIR
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpec
import wdlTools.generators.code.Fixer

import java.net.URI
import java.nio.file.{Files, Path, Paths}

class CorporaTest extends AnyWordSpec with Matchers {
  private val corporaDir = Option(getClass.getResource(s"/corpora")).map(f => Paths.get(f.getPath))
  private val configFile =
    Option(getClass.getResource(s"/corpora_repos.json")).map(f => Paths.get(f.getPath))
  private val corporaAvailable =
    configFile.exists(Files.exists(_)) && corporaDir.exists(Files.isDirectory(_))

  if (corporaAvailable) {
    val dxApi = DxApi()(Logger.Quiet)
    val dxProjectId = dxApi.currentProjectId.get
    val compilerOpts =
      Vector("--compileMode", "ir", "-quiet", "--locked", "--project", dxProjectId)
    lazy val fixDir = {
      val fixDir = Files.createTempDirectory("fix")
      fixDir.toFile.deleteOnExit()
      fixDir
    }
    lazy val uniqueDir = Iterator.from(0)

    def compile(dir: Path, file: Path, fix: Boolean, importDirs: Vector[Path]): Unit = {
      val path = if (fix) {
        // If fix is true, we need to fix the document first and write it to a tempdir.
        // The Fixer takes care of checking that the fixed file is valid.
        val fileResolver = FileSourceResolver.create(dir +: importDirs)
        val sourceFile = fileResolver.fromPath(dir.resolve(file))
        val fixer = Fixer(fileResolver = fileResolver)
        fixer.fix(sourceFile, fixDir.resolve(uniqueDir.next().toString)) match {
          case fs: LocalFileSource => fs.canonicalPath
          case _                   => throw new Exception("expected local file")
        }
      } else {
        dir.resolve(file)
      }
      Main.compile(path.toString +: compilerOpts) match {
        case SuccessfulCompileIR(_) => ()
        case err: UnsuccessfulTermination if err.exception.isDefined =>
          throw new Exception(s"failed to compile ${path}: ${err.message}", err.exception.get)
        case err: UnsuccessfulTermination =>
          throw new Exception(s"failed to compile ${path}: ${err.message}")
        case _ => throw new Exception("unexpected")
      }
    }

    val corpora = JsUtils.getValues(JsUtils.jsFromFile(configFile.get), Some("corpora"))

    "corpora test" should {
      corpora.map(JsUtils.getFields(_)).foreach { corpus =>
        val uri = URI.create(JsUtils.getString(corpus("url")))
        val repo = Paths.get(uri.getPath).getFileName.toString match {
          case s if s.endsWith(".git") => s.dropRight(4)
          case s                       => s
        }
        val root = corporaDir.get.resolve(repo)
        JsUtils.getValues(corpus("entrypoints")).map(JsUtils.getFields(_)).foreach { example =>
          val fix = JsUtils.getOptionalBoolean(example, "fix").getOrElse(false)
          val expectFailure = JsUtils.getOptionalBoolean(example, "fail").getOrElse(false)
          if (!expectFailure) {
            val importDirs = JsUtils
              .getOptionalValues(example, "import_dirs")
              .map(_.map(d => root.resolve(JsUtils.getString(d))))
              .getOrElse(Vector.empty)
            if (example.contains("path")) {
              val path = Paths.get(JsUtils.getString(example("path")))
              s"compile ${root.resolve(path)}" in {
                compile(root, path, fix, importDirs)
              }
            } else {
              val dir = JsUtils.getOptionalString(example, "dir").map(root.resolve).getOrElse(root)
              val include = JsUtils
                .getOptionalValues(example, "include")
                .map(_.map(i => Paths.get(JsUtils.getString(i))).toSet)
              val exclude = JsUtils
                .getOptionalValues(example, "exclude")
                .map(_.map(e => Paths.get(JsUtils.getString(e))).toSet)
              Files
                .list(dir)
                .filter(!Files.isDirectory(_))
                .filter { path =>
                  path.getFileName.toString.endsWith(".wdl") &&
                  include.forall(_.exists(i => path.endsWith(i))) &&
                  !exclude.exists(_.exists(e => path.endsWith(e)))
                }
                .forEach { path =>
                  s"compile ${dir.resolve(path)}" in {
                    compile(dir, path, fix, importDirs)
                  }
                }
            }
          }
        }
      }
    }
  }
}
