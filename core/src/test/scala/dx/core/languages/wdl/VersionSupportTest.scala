package dx.core.languages.wdl

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.syntax.SyntaxException

import java.nio.file.{Path, Paths}

class VersionSupportTest extends AnyFlatSpec with Matchers {

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  "VersionSupport" should "throw a relevant exception if wdl source doesn't have version declaration" in {
    val thrown = the[SyntaxException] thrownBy
      VersionSupport.fromSourceFile(pathFromBasename("bugs", "apps_912_version_error_msg_fix.wdl"))
    thrown.getCause.getMessage should include(
        "Draft 2 version of WDL should not contain a formal `input` section"
    )
  }
}
