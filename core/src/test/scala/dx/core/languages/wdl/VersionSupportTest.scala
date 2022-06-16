package dx.core.languages.wdl

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.{Path, Paths}

class VersionSupportTest extends AnyFlatSpec with Matchers {

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  "VersionSupport" should "generate a document with default None for optional type in nested workflow" in {
    val (doc, _, versionSupport) =
      VersionSupport.fromSourceFile(
          pathFromBasename("bugs", "apps_1222_optional_default_none_outer.wdl")
      )
    val generatedDoc = versionSupport.generateDocument(doc)
    generatedDoc shouldBe ()
  }
}
