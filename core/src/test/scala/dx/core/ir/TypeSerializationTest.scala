package dx.core.ir

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsString, JsValue}

import scala.collection.immutable.SeqMap

class TypeSerializationTest extends AnyFlatSpec with Matchers {
  private val testCases: Vector[Type] = Vector(
      // Primitive types
      Type.TBoolean,
      Type.TInt,
      Type.TFloat,
      Type.TString,
      Type.TFile,
      // arrays
      Type.TArray(Type.TString),
      Type.TArray(Type.TFile, nonEmpty = true),
      // optionals
      Type.TOptional(Type.TInt),
      Type.TOptional(Type.TArray(Type.TBoolean))
  )

  it should "work for various WDL types" in {
    testCases.foreach { t =>
      val typeMap = Map("key" -> t)
      val jsv = TypeSerde.serializeSpec(typeMap)
      TypeSerde.deserializeSpec(jsv, Map.empty) shouldBe typeMap
    }
  }

  private val personType =
    Type.TSchema("Person", SeqMap("name" -> Type.TString, "age" -> Type.TInt))
  private val houseType = Type.TSchema(
      "House",
      SeqMap("street" -> Type.TString, "zip code" -> Type.TInt, "owner" -> personType)
  )
  private val structTestCases: Vector[Type] = Vector(
      personType,
      Type.TArray(houseType),
      Type.TOptional(houseType)
  )

  it should "work for structs" in {
    val typeAliases: Map[String, Type.TSchema] = Map("Person" -> personType, "House" -> houseType)
    structTestCases.foreach { t =>
      val typeMap = Map("key" -> t)
      val jsv = TypeSerde.serializeSpec(typeMap)
      TypeSerde.deserializeSpec(jsv, typeAliases) shouldBe typeMap
    }
  }

  val badTypes: Vector[JsValue] = Vector(
      JsString("A bad type"),
      JsString("placeholder"),
      JsObject("type" -> JsString("Map"),
               "keys" -> JsString("Int"),
               "values" -> JsString("UnrealFile"))
  )

  it should "detect bad type descriptions" in {
    val typeAliases: Map[String, Type.TSchema] = Map("Person" -> personType, "House" -> houseType)
    badTypes.foreach { jsv =>
      assertThrows[Exception] {
        val typeMap = JsObject("key" -> jsv)
        TypeSerde.deserializeSpec(typeMap, typeAliases)
      }
    }
  }
}
