package dx.core.languages.wdl

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._
import wdlTools.eval.WdlValues._
import wdlTools.types.WdlTypes._

import scala.collection.immutable.TreeSeqMap

class WdlManifestTest extends AnyFlatSpec with Matchers {
  val simpleManifestJs: JsValue =
    """
      |{
      |  "person": {
      |    "name": {
      |      "first": "John",
      |      "last": "Smith"
      |    },
      |    "age": 42
      |  },
      |  "height_ft_in": [6, 1],
      |  "bank_account": "1234567890"
      |}""".stripMargin.parseJson

  val extendedManifestJs: JsValue =
    """
      |{
      |  "schemas": {
      |    "Name": {
      |      "type": "Struct",
      |      "members": {
      |        "first": "String",
      |        "last": "String"
      |      }
      |    },
      |    "Person": {
      |      "type": "Struct",
      |      "members": {
      |        "name": "Name",
      |        "age": "Int"
      |      }
      |    },
      |    "OptionalIntArray": {
      |      "type": "Array",
      |      "items": "Int",
      |      "nonEmpty": false,
      |      "optional": true
      |    }
      |  },
      |  "types": {
      |    "person": "Person",
      |    "bank_account": "String",
      |    "height_ft_in": {
      |      "type": "Pair",
      |      "left": "Int",
      |      "right": "Int"
      |    },
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
      |    "height_ft_in": [6, 1],
      |    "bank_account": "1234567890"
      |  }
      |}""".stripMargin.parseJson

  val serializedExtendedManifestJs: JsValue =
    """
      |{
      |  "schemas": {
      |    "Name": {
      |      "type": "Struct",
      |      "members": {
      |        "first": "String",
      |        "last": "String"
      |      }
      |    },
      |    "Person": {
      |      "type": "Struct",
      |      "members": {
      |        "age": "Int",
      |        "name": "Name"
      |      }
      |    }
      |  },
      |  "types": {
      |    "bank_account": "String",
      |    "height_ft_in": {
      |      "left": "Int",
      |      "right": "Int",
      |      "type": "Pair"
      |    },
      |    "person": "Person",
      |    "ssn_digits": {
      |      "items": "Int",
      |      "nonEmpty": false,
      |      "optional": true,
      |      "type": "Array"
      |    }
      |  },
      |  "values": {
      |    "bank_account": "1234567890",
      |    "height_ft_in": {
      |      "left": 6,
      |      "right": 1
      |    },
      |    "person": {
      |      "age": 42,
      |      "name": {
      |        "first": "John",
      |        "last": "Smith"
      |      }
      |    }
      |  }
      |}""".stripMargin.parseJson

  it should "parse a simple manifest" in {
    val nameType = T_Struct("Name", TreeSeqMap("first" -> T_String, "last" -> T_String))
    val personType = T_Struct("Person",
                              TreeSeqMap(
                                  "name" -> nameType,
                                  "age" -> T_Int
                              ))
    val types: Map[String, T] = Map(
        "person" -> personType,
        "bank_account" -> T_String,
        "height_ft_in" -> T_Pair(T_Int, T_Int),
        "ssn_digits" -> T_Optional(T_Array(T_Int, nonEmpty = false))
    )
    val manifestParser = WdlManifestParser()
    val manifest = manifestParser.parse(simpleManifestJs, Some(types))
    manifest.values should contain theSameElementsAs
      Map(
          "person" -> V_Struct("Person",
                               Map(
                                   "name" -> V_Struct("Name",
                                                      members = Map("first" -> V_String("John"),
                                                                    "last" -> V_String("Smith"))),
                                   "age" -> V_Int(42)
                               )),
          "height_ft_in" -> V_Pair(V_Int(6), V_Int(1)),
          "bank_account" -> V_String("1234567890")
      )
    manifest.schemas should contain theSameElementsAs
      Map(
          "Name" -> nameType,
          "Person" -> personType
      )
  }

  it should "parse an extended manifest" in {
    val manifestParser = WdlManifestParser()
    val manifest = manifestParser.parse(extendedManifestJs)
    manifest.values should contain theSameElementsAs
      Map(
          "person" -> V_Struct("Person",
                               Map(
                                   "name" -> V_Struct("Name",
                                                      members = Map("first" -> V_String("John"),
                                                                    "last" -> V_String("Smith"))),
                                   "age" -> V_Int(42)
                               )),
          "height_ft_in" -> V_Pair(V_Int(6), V_Int(1)),
          "bank_account" -> V_String("1234567890")
      )
    val nameType = T_Struct("Name", TreeSeqMap("first" -> T_String, "last" -> T_String))
    val personType = T_Struct("Person",
                              TreeSeqMap(
                                  "name" -> nameType,
                                  "age" -> T_Int
                              ))
    manifest.types should contain theSameElementsAs
      Map(
          "person" -> personType,
          "bank_account" -> T_String,
          "height_ft_in" -> T_Pair(T_Int, T_Int),
          "ssn_digits" -> T_Optional(T_Array(T_Int, nonEmpty = false))
      )
    manifest.schemas should contain theSameElementsAs
      Map(
          "Name" -> nameType,
          "Person" -> personType
      )
  }

  it should "serialize a manifest" in {
    val nameType = T_Struct("Name", TreeSeqMap("first" -> T_String, "last" -> T_String))
    val personType = T_Struct("Person",
                              TreeSeqMap(
                                  "name" -> nameType,
                                  "age" -> T_Int
                              ))
    val manifest = WdlManifest(
        TreeSeqMap(
            "person" -> V_Struct("Person",
                                 Map(
                                     "name" -> V_Struct("Name",
                                                        members = Map("first" -> V_String("John"),
                                                                      "last" -> V_String("Smith"))),
                                     "age" -> V_Int(42)
                                 )),
            "height_ft_in" -> V_Pair(V_Int(6), V_Int(1)),
            "bank_account" -> V_String("1234567890")
        ),
        TreeSeqMap(
            "person" -> personType,
            "bank_account" -> T_String,
            "height_ft_in" -> T_Pair(T_Int, T_Int),
            "ssn_digits" -> T_Optional(T_Array(T_Int, nonEmpty = false))
        ),
        TreeSeqMap(
            "Name" -> nameType,
            "Person" -> personType
        )
    )
    manifest.toJson.prettyPrint shouldBe serializedExtendedManifestJs.prettyPrint
  }
}
