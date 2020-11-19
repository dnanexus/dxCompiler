package dx.core.languages.cwl

import dx.cwl.{CwlFloat, CwlNull, CwlString, CwlType, FloatValue, StringValue}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class UtilsTest extends AnyFlatSpec with Matchers {
  it should "get default value given cwltype" in {
    val stringType = Vector[CwlType](CwlString)
    val floatType = Vector[CwlType](CwlFloat)
    val empty = Vector.empty[CwlType]
    val optionalStringType = Vector[CwlType](CwlNull, CwlString)

    CwlUtils.getDefaultCWLValue(stringType).shouldBe(StringValue(""))
    CwlUtils.getDefaultCWLValue(floatType).shouldBe(FloatValue(0.0))
    CwlUtils.getDefaultCWLValue(optionalStringType).shouldBe(StringValue(""))
    assertThrows[Exception](CwlUtils.getDefaultCWLValue(empty))
  }
}
