package dx.core.languages.wdl

import dx.Tags.ApiTest
import dx.api.{DiskType, DxInstanceType, InstanceTypeDB, InstanceTypeRequest}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._
import wdlTools.eval.{DefaultEvalPaths, Eval, Runtime => WdlRuntime}
import wdlTools.syntax.{SourceLocation, WdlVersion}
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import dx.util.{CodecUtils, FileSourceResolver, Logger}

import scala.collection.immutable.TreeSeqMap

class InstanceTypesTest extends AnyFlatSpec with Matchers {
  // The original list is at:
  // https://github.com/dnanexus/nucleus/blob/master/commons/instance_types/aws_instance_types.json
  //
  // The g2,i2,x1 instances have been removed, because they are not
  // enabled for customers by default.  In addition, the PV (Paravirtual)
  // instances have been removed, because they work only on Ubuntu
  // 12.04.
  //
  // Removed the ssd2 instances, because they actually use EBS storage. A better
  // solution would be asking the platform for the available instances.
  private val instanceList: String =
    """{
        "mem2_ssd1_x2":       {"internalName": "m3.large",                          "traits": {"numCores":   2, "totalMemoryMB":    7225, "ephemeralStorageGB":   27}},
        "mem2_ssd1_x4":       {"internalName": "m3.xlarge",                         "traits": {"numCores":   4, "totalMemoryMB":   14785, "ephemeralStorageGB":   72}},
        "mem2_ssd1_x8":       {"internalName": "m3.2xlarge",                        "traits": {"numCores":   8, "totalMemoryMB":   29905, "ephemeralStorageGB":  147}},
        "mem1_ssd1_x2":       {"internalName": "c3.large",                          "traits": {"numCores":   2, "totalMemoryMB":    3766, "ephemeralStorageGB":   28}},
        "mem1_ssd1_x4":       {"internalName": "c3.xlarge",                         "traits": {"numCores":   4, "totalMemoryMB":    7225, "ephemeralStorageGB":   77}},
        "mem1_ssd1_x8":       {"internalName": "c3.2xlarge",                        "traits": {"numCores":   8, "totalMemoryMB":   14785, "ephemeralStorageGB":  157}},
        "mem1_ssd1_x16":      {"internalName": "c3.4xlarge",                        "traits": {"numCores":  16, "totalMemoryMB":   29900, "ephemeralStorageGB":  302}},
        "mem1_ssd1_x32":      {"internalName": "c3.8xlarge",                        "traits": {"numCores":  32, "totalMemoryMB":   60139, "ephemeralStorageGB":  637}},
        "mem3_ssd1_x32_gen1": {"internalName": "cr1.8xlarge",                       "traits": {"numCores":  32, "totalMemoryMB":  245751, "ephemeralStorageGB":  237}},
        "mem3_ssd1_x2":       {"internalName": "r3.large",                          "traits": {"numCores":   2, "totalMemoryMB":   15044, "ephemeralStorageGB":   27}},
        "mem3_ssd1_x4":       {"internalName": "r3.xlarge",                         "traits": {"numCores":   4, "totalMemoryMB":   30425, "ephemeralStorageGB":   72}},
        "mem3_ssd1_x8":       {"internalName": "r3.2xlarge",                        "traits": {"numCores":   8, "totalMemoryMB":   61187, "ephemeralStorageGB":  147}},
        "mem3_ssd1_x16":      {"internalName": "r3.4xlarge",                        "traits": {"numCores":  16, "totalMemoryMB":  122705, "ephemeralStorageGB":  297}},
        "mem3_ssd1_x32":      {"internalName": "r3.8xlarge",                        "traits": {"numCores":  32, "totalMemoryMB":  245751, "ephemeralStorageGB":  597}}
}"""

  private val awsOnDemandHourlyPrice =
    """|{
       | "cc2.8xlarge": 2.000,
       | "cg1.4xlarge": 2.100,
       | "m4.large": 0.108,
       | "m4.xlarge": 0.215,
       | "m4.2xlarge": 0.431,
       | "m4.4xlarge": 0.862,
       | "m4.10xlarge": 2.155,
       | "m4.16xlarge": 3.447,
       | "c4.large": 0.100,
       | "c4.xlarge": 0.199,
       | "c4.2xlarge": 0.398,
       | "c4.4xlarge": 0.796,
       | "c4.8xlarge": 1.591,
       | "p2.xlarge": 0.900,
       | "p2.8xlarge": 7.200,
       | "p2.16xlarge": 14.400,
       | "g2.2xlarge": 0.650,
       | "g2.8xlarge": 2.600,
       | "x1.16xlarge": 6.669,
       | "x1.32xlarge": 13.338,
       | "r4.large": 0.133,
       | "r4.xlarge": 0.266,
       | "r4.2xlarge": 0.532,
       | "r4.4xlarge": 1.064,
       | "r4.8xlarge": 2.128,
       | "r4.16xlarge": 4.256,
       | "r3.large": 0.166,
       | "r3.xlarge": 0.333,
       | "r3.2xlarge": 0.665,
       | "r3.4xlarge": 1.330,
       | "r3.8xlarge": 2.660,
       | "i2.xlarge": 0.853,
       | "i2.2xlarge": 1.705,
       | "i2.4xlarge": 3.410,
       | "i2.8xlarge": 6.820,
       | "d2.xlarge": 0.690,
       | "d2.2xlarge": 1.380,
       | "d2.4xlarge": 2.760,
       | "d2.8xlarge": 5.520,
       | "hi1.4xlarge": 3.100,
       | "hs1.8xlarge": 4.600,
       | "m3.medium": 0.067,
       | "m3.large": 0.133,
       | "m3.xlarge": 0.266,
       | "m3.2xlarge": 0.532,
       | "c3.large": 0.090,
       | "c3.xlarge": 0.210,
       | "c3.2xlarge": 0.420,
       | "c3.4xlarge": 0.840,
       | "c3.8xlarge": 1.680,
       | "m1.small": 0.044,
       | "m1.medium": 0.087,
       | "m1.large": 0.175,
       | "m1.xlarge": 0.350,
       | "c1.medium": 0.130,
       | "c1.xlarge": 0.520,
       | "m2.xlarge": 0.245,
       | "m2.2xlarge": 0.490,
       | "m2.4xlarge": 0.980,
       | "t1.micro": 0.020,
       | "cr1.8xlarge": 3.500
       |}
       |""".stripMargin.trim

  // Create an availble instance list based on a hard coded list
  private def genTestDB(pricingInfo: Boolean): InstanceTypeDB = {
    def intOfJs(jsVal: JsValue): Int = {
      jsVal match {
        case JsNumber(x) => x.toInt
        case _           => throw new Exception("unexpected")
      }
    }
    val awsOnDemandHourlyPriceTable: Map[String, Float] = {
      val fields: Map[String, JsValue] = awsOnDemandHourlyPrice.parseJson.asJsObject.fields
      fields.map {
        case (name, v) =>
          val price: Float = v match {
            case JsNumber(x) => x.toFloat
            case _           => throw new Exception("unexpected")
          }
          name -> price
      }
    }

    val allInstances: Map[String, JsValue] = instanceList.parseJson.asJsObject.fields
    val db = allInstances.map {
      case (name, v) =>
        val fields: Map[String, JsValue] = v.asJsObject.fields
        val internalName = fields("internalName") match {
          case JsString(s) => s
          case _           => throw new Exception("unexpected")
        }
        val price: Option[Float] =
          if (pricingInfo) Some(awsOnDemandHourlyPriceTable(internalName))
          else None
        val traits = fields("traits").asJsObject.fields
        val memoryMB = intOfJs(traits("totalMemoryMB"))
        val diskGB = intOfJs(traits("ephemeralStorageGB"))
        val diskType = name.toLowerCase match {
          case s if s.contains("_ssd") => Some(DiskType.SSD)
          case s if s.contains("_hdd") => Some(DiskType.HDD)
          case _                       => None
        }
        val cpu = intOfJs(traits("numCores"))
        name -> DxInstanceType(name,
                               memoryMB,
                               diskGB,
                               cpu,
                               gpu = false,
                               Vector.empty,
                               diskType,
                               price)
    }
    InstanceTypeDB(db, pricingInfo)
  }

  private val dbFull = genTestDB(true)
  private lazy val dbOpaque = InstanceTypeDB.opaquePrices(dbFull)
  private val evaluator: Eval =
    Eval(DefaultEvalPaths.empty,
         Some(WdlVersion.V1),
         Vector.empty,
         FileSourceResolver.get,
         Logger.get)

  private def createRuntime(dxInstanceType: Option[String],
                            memory: Option[String],
                            disks: Option[String],
                            cpu: Option[String],
                            gpu: Option[Boolean]): Runtime = {
    def makeString(s: String): TAT.Expr = TAT.ValueString(s, WdlTypes.T_String, null)
    val rt = TreeSeqMap(
        Runtime.DxInstanceTypeKey -> dxInstanceType.map(makeString),
        WdlRuntime.MemoryKey -> memory.map(makeString),
        WdlRuntime.DisksKey -> disks.map(makeString),
        WdlRuntime.CpuKey -> cpu.map(makeString),
        WdlRuntime.GpuKey -> gpu.map(b => TAT.ValueBoolean(b, WdlTypes.T_Boolean, null))
    ).collect {
      case (key, Some(value)) => key -> value
    }
    Runtime(
        WdlVersion.V1,
        Some(TAT.RuntimeSection(rt, null)),
        None,
        evaluator
    )
  }

  private def useDB(db: InstanceTypeDB): Unit = {
    db.selectOptimal(InstanceTypeRequest.empty) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem1_ssd1_x2" =>
    }
    db.selectOptimal(
        InstanceTypeRequest(minMemoryMB = Some(3 * 1024), minDiskGB = Some(100), minCpu = Some(5))
    ) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem1_ssd1_x8" =>
    }
    db.selectOptimal(InstanceTypeRequest(minMemoryMB = Some(2 * 1024), minDiskGB = Some(20))) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem1_ssd1_x2" =>
    }
    db.selectOptimal(
        InstanceTypeRequest(minMemoryMB = Some(30 * 1024), minDiskGB = Some(128), minCpu = Some(8))
    ) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem3_ssd1_x8" =>
    }

    // no instance with 1024 CPUs
    db.selectOptimal(InstanceTypeRequest(minCpu = Some(1024))) shouldBe None

    db.apply(
          createRuntime(None, Some("3 GB"), Some("local-disk 10 HDD"), Some("1"), None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x2")
    db.apply(
          createRuntime(None, Some("37 GB"), Some("local-disk 10 HDD"), Some("6"), None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x8")
    db.apply(
          createRuntime(None, Some("2 GB"), Some("local-disk 100 HDD"), None, None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x8")
    db.apply(
          createRuntime(None, Some("2.1GB"), Some("local-disk 100 HDD"), None, None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x8")

    db.apply(
          createRuntime(Some("mem3_ssd1_x8"), None, None, None, None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x8")

    db.apply(
          createRuntime(None, Some("235 GB"), Some("local-disk 550 HDD"), Some("32"), None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x32")
    db.apply(
          createRuntime(Some("mem3_ssd1_x32"), None, None, None, None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x32")

    db.apply(
          createRuntime(None, None, None, Some("8"), None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x8")
  }

  it should "Choose reasonable platform instance types" in {
    useDB(dbFull)
  }

  it should "work even with opaque prices" taggedAs ApiTest in {
    useDB(dbOpaque)
  }

  it should "catch parsing errors" in {
    // TODO: due to coercion, an int value is able to be coerced to
    //  a string, so this doesn't actually throw an exception at
    //  parse time, though it will throw an exception when the value
    //  is not found in the instance type DB
//    assertThrows[Exception] {
//      // illegal request format
//      Runtime(
//          WdlVersion.V1,
//          Some(
//              TAT.RuntimeSection(
//                  Map(Runtime.DxInstanceTypeKey -> TAT.ValueInt(4, WdlTypes.T_Int, null)),
//                  null
//              )
//          ),
//          None,
//          evaluator,
//          WdlValueBindings.empty
//      ).parseInstanceType
//    }

    // note that wdlTools applies default values to memory, disk, and cpu

    // memory specification
    createRuntime(None, Some("230MB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(
          minMemoryMB = Some(Math.ceil((230d * 1000d * 1000d) / (1024d * 1024d)).toLong),
          minDiskGB = Some(1),
          minCpu = Some(1)
      )

    createRuntime(None, Some("230MiB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(230), minDiskGB = Some(1), minCpu = Some(1))

    createRuntime(None, Some("230GB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(
          minMemoryMB = Some(
              Math.ceil((230d * 1000d * 1000d * 1000d) / (1024d * 1024d)).toLong
          ),
          minDiskGB = Some(1),
          minCpu = Some(1)
      )

    createRuntime(None, Some("230GiB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(230 * 1024), minDiskGB = Some(1), minCpu = Some(1))

    createRuntime(None, Some("1000 TB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(
          minMemoryMB =
            Some(Math.ceil((1000d * 1000d * 1000d * 1000d * 1000d) / (1024d * 1024d)).toLong),
          minDiskGB = Some(1),
          minCpu = Some(1)
      )

    createRuntime(None, Some("1000 TiB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(1000L * 1024L * 1024L),
                          minDiskGB = Some(1),
                          minCpu = Some(1))

    assertThrows[Exception] {
      createRuntime(None, Some("230 44 34 GB"), None, None, None).parseInstanceType
    }
    assertThrows[Exception] {
      createRuntime(None, Some("230.x GB"), None, None, None).parseInstanceType
    }
    assertThrows[Exception] {
      createRuntime(None, Some("230.x GB"), None, None, None).parseInstanceType
    }
    // this one is fine - 230.3 is just rounded up to 231
//    assertThrows[Exception] {
//      createRuntime(None, Some("230.3"), None, None, None).parseInstanceType
//    }
    assertThrows[Exception] {
      createRuntime(None, Some("230 XXB"), None, None, None).parseInstanceType
    }

    // disk spec
    assertThrows[Exception] {
      createRuntime(None, None, Some("just give me a disk"), None, None).parseInstanceType
    }
    assertThrows[Exception] {
      createRuntime(None, None, Some("local-disk xxxx"), None, None).parseInstanceType
    }
//    assertThrows[Exception] {
//      createRuntime(None, None, Some(1024), None, None).parseInstanceType.get
//    }

    // cpu
    assertThrows[Exception] {
      createRuntime(None, None, None, Some("xxyy"), None).parseInstanceType
    }
    createRuntime(None, None, None, Some("1"), None).parseInstanceType shouldBe InstanceTypeRequest(
        minMemoryMB = Some(2048),
        minDiskGB = Some(1),
        minCpu = Some(1)
    )
    createRuntime(None, None, None, Some("1.2"), None).parseInstanceType shouldBe InstanceTypeRequest(
        minMemoryMB = Some(2048),
        minDiskGB = Some(1),
        minCpu = Some(2)
    )

//    assertThrows[Exception] {
//      createRuntime(None, None, None, Some(V_Boolean(false)), None)
//    }

    // gpu
    createRuntime(None, Some("1000 TiB"), None, None, Some(true)).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(1000L * 1024L * 1024L),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          gpu = Some(true))

    createRuntime(None, None, None, None, Some(false)).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(2048),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          gpu = Some(false))
  }

  it should "get required instance type" in {
    val encodedDb =
      "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"
    val js = CodecUtils.base64DecodeAndGunzip(encodedDb)
    val db = js.parseJson.convertTo[InstanceTypeDB]
    // query: memory=Some(2048) disk=Some(1) diskType=None cores=Some(5) gpu=None  os=None instancetype=None
    val runtime = Runtime(
        WdlVersion.V1,
        Some(
            TAT.RuntimeSection(
                TreeSeqMap("cpu" -> TAT.ValueInt(5, WdlTypes.T_Int, SourceLocation.empty)),
                SourceLocation.empty
            )
        ),
        None,
        evaluator,
        None,
        None
    )
    val request = runtime.parseInstanceType
    db.apply(request).name
  }
}
