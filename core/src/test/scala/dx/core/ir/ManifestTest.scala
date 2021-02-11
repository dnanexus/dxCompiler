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
    val types: Map[String, Type] = Map(
        "person" -> personType,
        "bank_account" -> TString,
        "height_ft_in" -> pairType,
        "ssn_digits" -> TOptional(TArray(TInt))
    )
    val manifestParser = ManifestParser()
    val manifest = manifestParser.parseUserManifest(userManifestJs, types)
    manifest.values should contain theSameElementsAs
      Map(
          "person" -> VHash(
              TreeSeqMap(
                  "name" -> VHash(
                      TreeSeqMap("first" -> VString("John"), "last" -> VString("Smith"))
                  ),
                  "age" -> VInt(42)
              )
          ),
          "height_ft_in" -> VHash(TreeSeqMap("left" -> VInt(6), "right" -> VInt(1))),
          "bank_account" -> VString("1234567890")
      )
    manifest.definitions should contain theSameElementsAs
      Map(
          "Name" -> nameType,
          "Person" -> personType,
          "Pair___(Int, Int)" -> pairType
      )
  }

  it should "parse a full manifest" in {
    val manifestParser = ManifestParser()
    val manifest = manifestParser.parseFullManifest(fullManifestJs)
    manifest.id shouldBe Some("test")
    manifest.values should contain theSameElementsAs
      Map(
          "person" -> VHash(
              TreeSeqMap(
                  "name" -> VHash(
                      TreeSeqMap("first" -> VString("John"), "last" -> VString("Smith"))
                  ),
                  "age" -> VInt(42)
              )
          ),
          "height_ft_in" -> VHash(TreeSeqMap("left" -> VInt(6), "right" -> VInt(1))),
          "bank_account" -> VString("1234567890")
      )
    val nameType = TSchema("Name", TreeSeqMap("first" -> TString, "last" -> TString))
    val personType = TSchema("Person", TreeSeqMap("name" -> nameType, "age" -> TInt))
    val pairType = TSchema("Pair___(Int, Int)", TreeSeqMap("left" -> TInt, "right" -> TInt))
    manifest.types should contain theSameElementsAs
      Map(
          "person" -> personType,
          "bank_account" -> TString,
          "height_ft_in" -> pairType,
          "ssn_digits" -> TOptional(TArray(TInt))
      )
    manifest.definitions should contain theSameElementsAs
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
    val manifest = TypedManifest(
        TreeSeqMap(
            "bank_account" -> VString("1234567890"),
            "height_ft_in" -> VHash(TreeSeqMap("left" -> VInt(6), "right" -> VInt(1))),
            "person" -> VHash(
                TreeSeqMap(
                    "age" -> VInt(42),
                    "name" -> VHash(
                        TreeSeqMap("first" -> VString("John"), "last" -> VString("Smith"))
                    )
                )
            )
        ),
        TreeSeqMap(
            "bank_account" -> TString,
            "height_ft_in" -> pairType,
            "person" -> personType,
            "ssn_digits" -> TOptional(TArray(TInt))
        ),
        TreeSeqMap(
            "Name" -> nameType,
            "Pair___(Int, Int)" -> pairType,
            "Person" -> personType
        ),
        id = Some("test")
    )
    manifest.toJson.prettyPrint shouldBe serializedManifestJs.prettyPrint
  }
}
