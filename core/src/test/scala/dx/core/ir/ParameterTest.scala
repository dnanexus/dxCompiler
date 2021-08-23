package dx.core.ir

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class ParameterTest extends AnyFlatSpec with Matchers {
  it should "not encode names with illegal characters" in {
    assertThrows[Exception] {
      Parameter.encodeName("  ")
    }
    assertThrows[Exception] {
      Parameter.encodeName("a__b")
    }
    assertThrows[Exception] {
      Parameter.encodeName("a b")
    }
    assertThrows[Exception] {
      Parameter.encodeName("foo.")
    }
    assertThrows[Exception] {
      Parameter.encodeName("a_._b")
    }
  }

  it should "correctly encode and decode names" in {
    Parameter.encodeName("a.b.c") shouldBe "a___b___c"
    Parameter.decodeName("a___b___c") shouldBe "a.b.c"

    Parameter.encodeName("_foo/bar_") shouldBe "_foo__47__bar_"
    Parameter.decodeName("_foo__47__bar_") shouldBe "_foo/bar_"

    Parameter.encodeName("a-9/c.d") shouldBe "a__45__9__47__c___d"
    Parameter.decodeName("a__45__9__47__c___d") shouldBe "a-9/c.d"

    Parameter.encodeName("foo.bar___") shouldBe "foo___bar___"
    Parameter.decodeName("foo___bar___") shouldBe "foo.bar___"
  }

  it should "correctly encode and decode simple names" in {
    Parameter.encodeDots("a.b.c") shouldBe "a___b___c"
    Parameter.decodeDots("a___b___c") shouldBe "a.b.c"

    Parameter.encodeDots("_foo.bar_") shouldBe "_foo___bar_"
    Parameter.decodeDots("_foo___bar_") shouldBe "_foo.bar_"

    Parameter.encodeDots("foo.bar___") shouldBe "foo___bar___"
    Parameter.decodeDots("foo___bar___") shouldBe "foo.bar___"
  }
}
