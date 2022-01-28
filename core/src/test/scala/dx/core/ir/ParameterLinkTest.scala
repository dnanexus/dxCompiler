package dx.core.ir

import java.nio.file.Paths
import dx.core.Tags.EdgeTest
import dx.api.{DxApi, DxFileDescCache}
import dx.core.Constants
import dx.core.languages.wdl.WdlDxName
import dx.util.{FileSourceResolver, Logger}
import dx.util.protocols.DxFileAccessProtocol
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

import scala.collection.immutable.SeqMap

class ParameterLinkTest extends AnyFlatSpec with Matchers {
  private val dxApi = DxApi()(Logger.Quiet)
  private val dxProtocol = DxFileAccessProtocol(dxApi)
  private val fileResolver = FileSourceResolver.create(userProtocols = Vector(dxProtocol))
  private val parameterLinkSerializer = ParameterLinkSerializer(fileResolver, dxApi)
  private val parameterLinkDeserializer = ParameterLinkDeserializer(DxFileDescCache.empty, dxApi)

  case class Element(name: String, irType: Type, irValue: Value)

  def makeElement(t: Type, v: Value): Element = Element("A", t, v)

  def check(elem: Element): Unit = {
    val prefix = "XXX_"
    val link = parameterLinkSerializer.createLink(elem.irType, elem.irValue)
    val dxName = WdlDxName.fromSourceName(prefix + elem.name)
    val allDxFields1: Vector[(DxName, JsValue)] =
      parameterLinkSerializer.createFieldsFromLink(link, dxName)
    val allDxFields2 = allDxFields1.filter {
      case (key, _) => !key.suffix.contains(Constants.FlatFilesSuffix)
    }
    allDxFields2.size should be(1)
    val (name2, jsv) = allDxFields2.head
    name2 shouldBe dxName
    val wdlValue2 =
      parameterLinkDeserializer.deserializeInputWithType(jsv, elem.irType, dxName.decoded)
    wdlValue2 should be(elem.irValue)
  }

  it should "handle primitive WDL elements" in {
    val testCases = Vector(
        // primitives
        makeElement(Type.TBoolean, Value.VBoolean(true)),
        makeElement(Type.TInt, Value.VInt(19)),
        makeElement(Type.TFloat, Value.VFloat(2.718)),
        makeElement(Type.TString, Value.VString("water and ice")),
        makeElement(Type.TFile, Value.VFile("/usr/var/local/bin/gcc"))
    )
    testCases.foreach(check)
  }

  it should "handle compound WDL types" in {
    val testCases = Vector(
        // optional
        makeElement(Type.TOptional(Type.TFile),
                    Value.VFile(Paths.get("ddd").toAbsolutePath.toString)),
        // arrays
        makeElement(
            Type.TArray(Type.TBoolean),
            Value.VArray(Vector(Value.VBoolean(true), Value.VBoolean(false)))
        ),
        // objects
        makeElement(
            Type.THash,
            Value.VHash(
                SeqMap(
                    "A" -> Value.VBoolean(true),
                    "C" -> Value.VBoolean(false),
                    "G" -> Value.VBoolean(true),
                    "H" -> Value.VBoolean(false)
                )
            )
        )
    )
    testCases.foreach(check)
  }

  it should "handle structs" in {
    val personType =
      Type.TSchema("Person", SeqMap("name" -> Type.TString, "age" -> Type.TInt))
    val jeff = Value.VHash(SeqMap("name" -> Value.VString("Jeoffrey"), "age" -> Value.VInt(16)))
    val janice = Value.VHash(SeqMap("name" -> Value.VString("Janice"), "age" -> Value.VInt(25)))
    val testCases = Vector(makeElement(personType, jeff), makeElement(personType, janice))

    // no definitions for struct Person, should fail
    // val wvlConverterEmpty = new WdlVarLinksConverter(verbose, Map.empty, Map.empty)
    // testCases.foreach{ elem =>
    // assertThrows[Exception] {
//                check(elem, wvlConverterEmpty)
//            }
//        }

    testCases.foreach(check)
  }

  it should "handle nested structs" taggedAs EdgeTest in {
    // People
    val personType =
      Type.TSchema("Person", SeqMap("name" -> Type.TString, "age" -> Type.TInt))
    val houseType = Type.TSchema(
        "House",
        SeqMap("person" -> personType, "zipcode" -> Type.TInt, "type" -> Type.TString)
    )

    // people
    val lucy = Value.VHash(SeqMap("name" -> Value.VString("Lucy"), "age" -> Value.VInt(37)))
    val lear =
      Value.VHash(SeqMap("name" -> Value.VString("King Lear"), "age" -> Value.VInt(41)))

    // Houses
    val learCastle =
      Value.VHash(
          SeqMap("person" -> lear, "zipcode" -> Value.VInt(1), "type" -> Value.VString("Castle"))
      )

    val lucyHouse =
      Value.VHash(
          SeqMap("person" -> lucy,
                 "zipcode" -> Value.VInt(94043),
                 "type" -> Value.VString("town house"))
      )

    val testCases = Vector(makeElement(houseType, learCastle), makeElement(houseType, lucyHouse))
    testCases.foreach(check)
  }
}
