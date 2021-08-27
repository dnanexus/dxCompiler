package dx.core.ir

import dx.core.ir.Type._
import dx.core.ir.Value._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

import scala.collection.immutable.TreeSeqMap

class ManifestTest extends AnyFlatSpec with Matchers {
  val userManifestJs: JsValue =
    """
      |{
      |  "person": {
      |    "name": {
      |      "first": "John",
      |      "last": "Smith"
      |    },
      |    "age": 42
      |  },
      |  "height_ft_in": {
      |    "left": 6, 
      |    "right": 1
      |  },
      |  "bank_account": "1234567890"
      |}""".stripMargin.parseJson

  val fullManifestJs: JsValue =
    """
      |{
      |  "id": "test",
      |  "definitions": {
      |    "Name": {
      |      "fields": {
      |        "first": "String",
      |        "last": "String"
      |      },
      |      "type": "Name"
      |    },
      |    "OptionalIntArray": {
      |      "items": "Int",
      |      "nonEmpty": false,
      |      "optional": true,
      |      "type": "Array"
      |    },
      |    "Pair___(Int, Int)": {
      |      "fields": {
      |        "left": "Int",
      |        "right": "Int"
      |      },
      |      "type": "Pair___(Int, Int)"
      |    },
      |    "Person": {
      |      "fields": {
      |        "name": "Name",
      |        "age": "Int"
      |      },
      |      "type": "Person"
      |    }
      |  },
      |  "types": {
      |    "person": "Person",
      |    "bank_account": "String",
      |    "height_ft_in": "Pair___(Int, Int)",
      |    "ssn_digits": "OptionalIntArray"
      |  },
      |  "values": {
      |    "person": {
      |      "name": {
      |        "first": "John",
      |        "last": "Smith"
      |      },
      |      "age": 42
      |    },
      |    "height_ft_in": {
      |      "left": 6, 
      |      "right": 1
      |    },
      |    "bank_account": "1234567890"
      |  }
      |}""".stripMargin.parseJson

  val serializedManifestJs: JsValue =
    """
      |{
      |  "definitions": {
      |    "Name": {
      |      "fields": {
      |        "first": "String",
      |        "last": "String"
      |      },
      |      "type": "Name"
      |    },
      |    "Pair___(Int, Int)": {
      |      "fields": {
      |        "left": "Int",
      |        "right": "Int"
      |      },
      |      "type": "Pair___(Int, Int)"
      |    },
      |    "Person": {
      |      "fields": {
      |        "name": "Name",
      |        "age": "Int"
      |      },
      |      "type": "Person"
      |    }
      |  },
      |  "encoded": false,
      |  "id": "test",
      |  "types": {
      |    "person": "Person",
      |    "bank_account": "String",
      |    "height_ft_in": "Pair___(Int, Int)",
      |    "ssn_digits": {
      |      "items": "Int",
      |      "nonEmpty": false,
      |      "optional": true,
      |      "type": "Array"
      |    }
      |  },
      |  "values": {
      |    "person": {
      |      "name": {
      |        "first": "John",
      |        "last": "Smith"
      |      },
      |      "age": 42
      |    },
      |    "height_ft_in": {
      |      "left": 6, 
      |      "right": 1
      |    },
      |    "bank_account": "1234567890"
      |  }
      |}""".stripMargin.parseJson

  it should "parse a user manifest" in {
    val nameType = TSchema("Name", TreeSeqMap("first" -> TString, "last" -> TString))
    val personType = TSchema("Person",
                             TreeSeqMap(
                                 "name" -> nameType,
                                 "age" -> TInt
                             ))
    val pairType = TSchema("Pair___(Int, Int)", TreeSeqMap("left" -> TInt, "right" -> TInt))
    val types: Map[DxName, Type] = Map(
        SimpleDxName.fromSourceName("person") -> personType,
        SimpleDxName.fromSourceName("bank_account") -> TString,
        SimpleDxName.fromSourceName("height_ft_in") -> pairType,
        SimpleDxName.fromSourceName("ssn_digits") -> TOptional(TArray(TInt))
    )
    val manifest = Manifest.parse(userManifestJs, SimpleDxName, Some(types))
    manifest.deserialize() should contain theSameElementsAs
      Map(
          SimpleDxName.fromSourceName("person") -> VHash(
              TreeSeqMap(
                  "name" -> VHash(
                      TreeSeqMap("first" -> VString("John"), "last" -> VString("Smith"))
                  ),
                  "age" -> VInt(42)
              )
          ),
          SimpleDxName.fromSourceName("height_ft_in") -> VHash(
              TreeSeqMap("left" -> VInt(6), "right" -> VInt(1))
          ),
          SimpleDxName.fromSourceName("bank_account") -> VString("1234567890")
      )
    manifest.definitions.nonEmpty shouldBe true
    manifest.definitions.get should contain theSameElementsAs
      Map(
          "Name" -> nameType,
          "Person" -> personType,
          "Pair___(Int, Int)" -> pairType
      )
  }

  it should "parse a full manifest" in {
    val manifest = Manifest.parse(fullManifestJs, SimpleDxName)
    manifest.id shouldBe Some("test")
    manifest.deserialize() should contain theSameElementsAs
      Map(
          SimpleDxName.fromSourceName("person") -> VHash(
              TreeSeqMap(
                  "name" -> VHash(
                      TreeSeqMap("first" -> VString("John"), "last" -> VString("Smith"))
                  ),
                  "age" -> VInt(42)
              )
          ),
          SimpleDxName.fromSourceName("height_ft_in") -> VHash(
              TreeSeqMap("left" -> VInt(6), "right" -> VInt(1))
          ),
          SimpleDxName.fromSourceName("bank_account") -> VString("1234567890")
      )
    val nameType = TSchema("Name", TreeSeqMap("first" -> TString, "last" -> TString))
    val personType = TSchema("Person", TreeSeqMap("name" -> nameType, "age" -> TInt))
    val pairType = TSchema("Pair___(Int, Int)", TreeSeqMap("left" -> TInt, "right" -> TInt))
    manifest.types.nonEmpty shouldBe true
    manifest.types.get should contain theSameElementsAs
      Map(
          SimpleDxName.fromSourceName("person") -> personType,
          SimpleDxName.fromSourceName("bank_account") -> TString,
          SimpleDxName.fromSourceName("height_ft_in") -> pairType,
          SimpleDxName.fromSourceName("ssn_digits") -> TOptional(TArray(TInt))
      )
    manifest.definitions.nonEmpty shouldBe true
    manifest.definitions.get should contain theSameElementsAs
      Map(
          "Name" -> nameType,
          "Person" -> personType,
          "Pair___(Int, Int)" -> pairType
      )
  }

  it should "serialize a manifest" in {
    val nameType = TSchema("Name", TreeSeqMap("first" -> TString, "last" -> TString))
    val personType = TSchema("Person",
                             TreeSeqMap(
                                 "name" -> nameType,
                                 "age" -> TInt
                             ))
    val pairType = TSchema("Pair___(Int, Int)", TreeSeqMap("left" -> TInt, "right" -> TInt))
    val manifest = Manifest(
        TreeSeqMap(
            SimpleDxName.fromSourceName("bank_account") -> JsString("1234567890"),
            SimpleDxName.fromSourceName("height_ft_in") -> JsObject("left" -> JsNumber(6),
                                                                    "right" -> JsNumber(1)),
            SimpleDxName.fromSourceName("person") -> JsObject(
                "age" -> JsNumber(42),
                "name" -> JsObject(
                    "first" -> JsString("John"),
                    "last" -> JsString("Smith")
                )
            )
        ),
        Some(
            TreeSeqMap(
                SimpleDxName.fromSourceName("bank_account") -> TString,
                SimpleDxName.fromSourceName("height_ft_in") -> pairType,
                SimpleDxName.fromSourceName("person") -> personType,
                SimpleDxName.fromSourceName("ssn_digits") -> TOptional(TArray(TInt))
            )
        ),
        Some(
            TreeSeqMap(
                "Name" -> nameType,
                "Pair___(Int, Int)" -> pairType,
                "Person" -> personType
            )
        ),
        id = Some("test")
    )
    manifest.toJson().prettyPrint shouldBe serializedManifestJs.prettyPrint
  }
}
