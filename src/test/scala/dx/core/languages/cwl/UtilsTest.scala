package dx.core.languages.cwl

import dx.core.ir.Type
import dx.core.ir.Type.{TFloat, TOptional, TString}
import dx.core.ir.Value.{VFloat, VInt, VString}
import dx.cwl.{CwlFloat, CwlNull, CwlString, CwlType, CwlValue, FloatValue, IntValue, StringValue}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class UtilsTest extends AnyFlatSpec with Matchers {
  val stringType: Vector[CwlType] = Vector[CwlType](CwlString)
  val floatType: Vector[CwlType] = Vector[CwlType](CwlFloat)
  val emptyType: Vector[CwlType] = Vector.empty[CwlType]
  val optionalStringType: Vector[CwlType] = Vector[CwlType](CwlNull, CwlString)
  val stringIRType: Type = TString
  val floatIRType: Type = TFloat
  val optionalIRType: Type = TOptional(TString)
  val stringCwl: CwlValue = StringValue("test_string_value")
  val intCwl: CwlValue = IntValue(3)

  it should "get default value given cwltype" in {
    CwlUtils.getDefaultCWLValue(stringType).shouldBe(StringValue(""))
    CwlUtils.getDefaultCWLValue(floatType).shouldBe(FloatValue(0.0))
    CwlUtils.getDefaultCWLValue(optionalStringType).shouldBe(StringValue(""))
    assertThrows[Exception](CwlUtils.getDefaultCWLValue(emptyType))
  }

  it should "translate CWLTypes to IRType" in {
    CwlUtils.toIRType(stringType).shouldBe(TString)
    CwlUtils.toIRType(floatType).shouldBe(TFloat)
    CwlUtils.toIRType(optionalStringType).shouldBe(TOptional(TString))
    assertThrows[Exception](CwlUtils.toIRType(emptyType))
  }

  it should "get default value from IRType" in {
    CwlUtils.getDefaultIRValue(stringType).shouldBe(VString(""))
    CwlUtils.getDefaultIRValue(stringIRType).shouldBe(VString(""))
    CwlUtils.getDefaultIRValue(floatType).shouldBe(VFloat(0.0))
    CwlUtils.getDefaultIRValue(floatIRType).shouldBe(VFloat(0.0))
    CwlUtils.getDefaultIRValue(optionalStringType).shouldBe(VString(""))
    CwlUtils.getDefaultIRValue(optionalIRType).shouldBe(VString(""))
  }

    it should "translate CwlValue to IRValue" in {
    CwlUtils.toIRValue(stringCwl).shouldBe(VString("test_string_value"))
    CwlUtils.toIRValue(intCwl).shouldBe(VInt(3))
  }
}
