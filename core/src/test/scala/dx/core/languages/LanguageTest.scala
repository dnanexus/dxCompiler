package dx.core.languages

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class LanguageTest extends AnyFlatSpec with Matchers {
  it should "parse various language strings" in {
    Language.parse("wdl 1.0") shouldBe Language.WdlV1_0
    Language.parse("wdl v2.0") shouldBe Language.WdlV2_0
    Language.parse("wdl") shouldBe Language.WdlDefault
    Language.parse("draft2", Some("wdl")) shouldBe Language.WdlVDraft2
    Language.parse("cwl") shouldBe Language.CwlDefault
    Language.parse("cwl v1.2") shouldBe Language.CwlV1_2
  }

  it should "fail to parse unsupported language strings" in {
    var thrown = intercept[Exception] {
      Language.parse("v1.0", Some("cwl"))
    }
    thrown.getMessage() should include (
      s"Unrecognized/unsupported language"
    )
    
    thrown = intercept[Exception] {
      Language.parse("v1.1", Some("cwl"))
    }
    thrown.getMessage() should include (
      s"Unrecognized/unsupported language"
    )
    
    thrown = intercept[Exception] {
      Language.parse("cwl1.0", Some("cwl"))
    }
    thrown.getMessage() should include (
      s"Unrecognized/unsupported language"
    )
  }
}
