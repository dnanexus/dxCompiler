package dx.core.ir

import dx.core.Constants
import dx.core.Constants.ComplexValueKey
import dx.core.languages.cwl.CwlDxName
import dx.core.languages.wdl.WdlDxName
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.util.Random

class DxNameTest extends AnyFlatSpec with Matchers {
  it should "split parts of a name" in {
    DxNameFactory.split("a.b.c______dxfiles", Some("\\.".r)) shouldBe (
        Vector("a", "b", "c"), None, Some("______dxfiles")
    )
  }

  it should "handle stage prefixes" in {
    DxNameFactory.split("stage-1.foo.bar___dxfiles", Some("\\.".r)) shouldBe (
        Vector("foo", "bar"), Some("stage-1"), Some("___dxfiles")
    )
  }

  it should "add a suffix" in {
    val dxName = SimpleDxName.fromSourceName("foo", Some(ComplexValueKey))
    dxName.encoded shouldBe "foo___"
    dxName.addSuffix(Constants.FlatFilesSuffix).encoded shouldBe "foo______dxfiles"
  }

  it should "encode and decode WDL names" in {
    WdlDxName.fromDecodedName("a.b.c").encoded shouldBe "a___b___c"
    WdlDxName.fromEncodedName("a___b___c").decoded shouldBe "a.b.c"

    WdlDxName.fromDecodedName("_foo.bar_").encoded shouldBe "_foo___bar_"
    WdlDxName.fromEncodedName("_foo___bar_").decoded shouldBe "_foo.bar_"

    WdlDxName.fromDecodedName("foo.bar___").encoded shouldBe "foo___bar___"
    WdlDxName.fromEncodedName("foo___bar___").decoded shouldBe "foo.bar___"

    WdlDxName.fromDecodedName("stage-1.foo.bar___").encoded shouldBe "stage-1.foo___bar___"
    WdlDxName.fromEncodedName("stage-1.foo___bar___").decoded shouldBe "stage-1.foo.bar___"
  }

  it should "encode a valid WDL name" in {
    val dxName = WdlDxName.fromDecodedName("stage-1.a.b___dxfiles")
    dxName.numParts shouldBe 2
    dxName.getDecodedParts shouldBe Vector("a", "b")
    dxName.getEncodedParts shouldBe Vector("a", "b")
    dxName.decoded shouldBe "stage-1.a.b___dxfiles"
    dxName.encoded shouldBe "stage-1.a___b___dxfiles"
    dxName.stage shouldBe Some("stage-1")
    dxName.suffix shouldBe Some("___dxfiles")

    val dxName2 = WdlDxName.fromDecodedName("b___dxfiles")
    dxName.endsWith(dxName2) shouldBe true

    val dxName3 = WdlDxName.fromDecodedName("stage-1.b___dxfiles")
    dxName3.pushDecodedNamespace("a") shouldBe dxName
    dxName.popDecodedIdentifier() shouldBe (WdlDxName.fromDecodedName("stage-1.a"), "b")
    dxName.popDecodedIdentifier(keepSuffix = true) shouldBe (
        WdlDxName.fromDecodedName("stage-1.a___dxfiles"), "b"
    )

    val dxName4 = dxName.dropSuffix()
    dxName4.suffix shouldBe None
    dxName4.encoded shouldBe "stage-1.a___b"
    dxName4.decoded shouldBe "stage-1.a.b"

    dxName4.withSuffix("___dxfiles") shouldBe dxName
  }

  it should "decode valid WDL names" in {
    val dxName = WdlDxName.fromEncodedName("a___b___dxfiles")
    dxName.numParts shouldBe 2
    dxName.getDecodedParts shouldBe Vector("a", "b")
    dxName.getEncodedParts shouldBe Vector("a", "b")
    dxName.decoded shouldBe "a.b___dxfiles"
    dxName.encoded shouldBe "a___b___dxfiles"
    dxName.suffix shouldBe Some("___dxfiles")

    val dxName2 = WdlDxName.fromEncodedName("b___dxfiles")
    dxName.endsWith(dxName2) shouldBe true
    dxName2.pushDecodedNamespace("a") shouldBe dxName
    dxName.popDecodedIdentifier() shouldBe (WdlDxName.fromEncodedName("a"), "b")

    val dxName3 = dxName.dropSuffix()
    dxName3.suffix shouldBe None
    dxName3.encoded shouldBe "a___b"
    dxName3.decoded shouldBe "a.b"

    dxName3.withSuffix("___dxfiles") shouldBe dxName
  }

  it should "not encode WDL names with illegal characters" in {
    assertThrows[Throwable] {
      WdlDxName.fromDecodedName("  ")
    }
    assertThrows[Throwable] {
      WdlDxName.fromDecodedName("a b")
    }
    assertThrows[Throwable] {
      WdlDxName.fromDecodedName("foo.")
    }
    assertThrows[Throwable] {
      WdlDxName.fromDecodedName("a_._b")
    }
    assertThrows[Throwable] {
      WdlDxName.fromDecodedName("a-b.c")
    }
  }

  it should "not allow encoded WDL names with illegal characters" in {
    assertThrows[Throwable] {
      WdlDxName.fromEncodedName("  ")
    }
    assertThrows[Throwable] {
      WdlDxName.fromEncodedName("a.")
    }
    assertThrows[Throwable] {
      WdlDxName.fromEncodedName("a.b")
    }
  }

  it should "sort WDL names correctly" in {
    val names = Vector(
        "a",
        "a.b",
        "a.b___dxfiles",
        "b"
    )
    (1 to 10).foreach { _ =>
      Random.shuffle(names).map(WdlDxName.fromDecodedName).sorted.map(_.decoded) shouldBe names
    }
  }

  it should "hash and compare two WDL names" in {
    val dxName1 = WdlDxName.fromEncodedName("a___b___dxfiles")
    val dxName2 = WdlDxName.fromDecodedName("a.b___dxfiles")
    (dxName1 == dxName2) shouldBe true
    val m = Map(dxName1 -> "x")
    m(dxName2) shouldBe "x"
  }

  it should "encode and decode CWL names" in {
    CwlDxName.fromDecodedName("a/b/c").encoded shouldBe "a___b___c"
    CwlDxName.fromEncodedName("a___b___c").decoded shouldBe "a/b/c"

    CwlDxName.fromDecodedName("_foo/bar_").encoded shouldBe "_foo___bar_"
    CwlDxName.fromEncodedName("_foo___bar_").decoded shouldBe "_foo/bar_"

    CwlDxName.fromDecodedName("a-9/c.d").encoded shouldBe "a__45__9___c__46__d"
    CwlDxName.fromEncodedName("a__45__9___c__46__d").decoded shouldBe "a-9/c.d"

    CwlDxName.fromDecodedName("foo/bar___").encoded shouldBe "foo___bar___"
    CwlDxName.fromEncodedName("foo___bar___").decoded shouldBe "foo/bar___"
  }

  it should "encode a valid CWL name" in {
    val dxName = CwlDxName.fromDecodedName("a-b/c.d___dxfiles")
    dxName.numParts shouldBe 2
    dxName.getDecodedParts shouldBe Vector("a-b", "c.d")
    dxName.getEncodedParts shouldBe Vector("a__45__b", "c__46__d")
    dxName.decoded shouldBe "a-b/c.d___dxfiles"
    dxName.encoded shouldBe "a__45__b___c__46__d___dxfiles"
    dxName.suffix shouldBe Some("___dxfiles")

    val dxName2 = CwlDxName.fromDecodedName("c.d___dxfiles")
    dxName.endsWith(dxName2) shouldBe true
    dxName2.pushDecodedNamespace("a-b") shouldBe dxName
    dxName.popDecodedIdentifier() shouldBe (CwlDxName.fromDecodedName("a-b"), "c.d")

    val dxName3 = dxName.dropSuffix()
    dxName3.suffix shouldBe None
    dxName3.encoded shouldBe "a__45__b___c__46__d"
    dxName3.decoded shouldBe "a-b/c.d"

    dxName3.withSuffix("___dxfiles") shouldBe dxName
  }

  it should "decode a valid CWL name" in {
    val dxName = CwlDxName.fromEncodedName("a__45__b___c__46__d___dxfiles")
    dxName.numParts shouldBe 2
    dxName.getDecodedParts shouldBe Vector("a-b", "c.d")
    dxName.getEncodedParts shouldBe Vector("a__45__b", "c__46__d")
    dxName.decoded shouldBe "a-b/c.d___dxfiles"
    dxName.encoded shouldBe "a__45__b___c__46__d___dxfiles"
    dxName.suffix shouldBe Some("___dxfiles")

    val dxName2 = CwlDxName.fromEncodedName("c__46__d___dxfiles")
    dxName.endsWith(dxName2) shouldBe true
    dxName2.pushDecodedNamespace("a-b") shouldBe dxName
    dxName.popDecodedIdentifier() shouldBe (CwlDxName.fromEncodedName("a__45__b"), "c.d")

    val dxName3 = dxName.dropSuffix()
    dxName3.suffix shouldBe None
    dxName3.encoded shouldBe "a__45__b___c__46__d"
    dxName3.decoded shouldBe "a-b/c.d"

    dxName3.withSuffix("___dxfiles") shouldBe dxName
  }

  it should "not allow decoded CWL names with illegal characters" in {
    assertThrows[Throwable] {
      CwlDxName.fromDecodedName("  ")
    }
    assertThrows[Throwable] {
      CwlDxName.fromDecodedName("a__b")
    }
    assertThrows[Throwable] {
      CwlDxName.fromDecodedName("a b")
    }
    assertThrows[Throwable] {
      CwlDxName.fromDecodedName("foo/")
    }
    assertThrows[Throwable] {
      CwlDxName.fromDecodedName("a_/_b")
    }
    assertThrows[Throwable] {
      CwlDxName.fromSourceName("a/b")
    }
  }

  it should "not allow encoded CWL names with illegal characters" in {
    assertThrows[Throwable] {
      CwlDxName.fromEncodedName("  ")
    }
    assertThrows[Throwable] {
      CwlDxName.fromEncodedName("a/")
    }
    assertThrows[Throwable] {
      CwlDxName.fromEncodedName("a/b")
    }
  }

  it should "sort CWL names correctly" in {
    val names = Vector(
        "a",
        "a-b/c.d",
        "a.b/c-d",
        "a/b",
        "a/b___dxfiles",
        "b"
    )
    (1 to 10).foreach { _ =>
      Random.shuffle(names).map(CwlDxName.fromDecodedName).sorted.map(_.decoded) shouldBe names
    }
  }

  it should "hash and compare two CWL names" in {
    val dxName1 = CwlDxName.fromEncodedName("a___b___dxfiles")
    val dxName2 = CwlDxName.fromDecodedName("a/b___dxfiles")
    (dxName1 == dxName2) shouldBe true
    val m = Map(dxName1 -> "x")
    m(dxName2) shouldBe "x"
  }
}
