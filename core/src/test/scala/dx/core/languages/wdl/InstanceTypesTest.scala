package dx.core.languages.wdl

import dx.api.{DiskType, DxApi, DxInstanceType, InstanceTypeDB, InstanceTypeRequest}

import dx.core.Constants
import dx.core.ir.RunSpec.InstanceType
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._
import wdlTools.eval.{DefaultEvalPaths, Eval, Runtime => WdlRuntime}
import wdlTools.syntax.{Quoting, SourceLocation, WdlVersion}
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}
import dx.util.{FileSourceResolver, JsUtils, Logger}

import scala.collection.immutable.SeqMap

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
        "mem2_ssd1_x2":       {"internalName": "m3.large",   "traits": {"numCores":   2, "totalminMemoryMB":    7225, "ephemeralStorageGB":   27, "rank": 2}},
        "mem2_ssd1_x4":       {"internalName": "m3.xlarge",  "traits": {"numCores":   4, "totalminMemoryMB":   14785, "ephemeralStorageGB":   72, "rank": 5}},
        "mem2_ssd1_x8":       {"internalName": "m3.2xlarge", "traits": {"numCores":   8, "totalminMemoryMB":   29905, "ephemeralStorageGB":  147, "rank": 8}},
        "mem1_ssd1_x2":       {"internalName": "c3.large",   "traits": {"numCores":   2, "totalminMemoryMB":    3766, "ephemeralStorageGB":   28, "rank": 1}},
        "mem1_ssd1_x4":       {"internalName": "c3.xlarge",  "traits": {"numCores":   4, "totalminMemoryMB":    7225, "ephemeralStorageGB":   77, "rank": 4}},
        "mem1_ssd1_x8":       {"internalName": "c3.2xlarge", "traits": {"numCores":   8, "totalminMemoryMB":   14785, "ephemeralStorageGB":  157, "rank": 7}},
        "mem1_ssd1_x16":      {"internalName": "c3.4xlarge", "traits": {"numCores":  16, "totalminMemoryMB":   29900, "ephemeralStorageGB":  302, "rank": 10}},
        "mem1_ssd1_x32":      {"internalName": "c3.8xlarge", "traits": {"numCores":  32, "totalminMemoryMB":   60139, "ephemeralStorageGB":  637, "rank": 12}},
        "mem3_ssd1_x32_gen1": {"internalName": "cr1.8xlarge","traits": {"numCores":  32, "totalminMemoryMB":  245751, "ephemeralStorageGB":  237, "rank": 14}},
        "mem3_ssd1_x2":       {"internalName": "r3.large",   "traits": {"numCores":   2, "totalminMemoryMB":   15044, "ephemeralStorageGB":   27, "rank": 3}},
        "mem3_ssd1_x4":       {"internalName": "r3.xlarge",  "traits": {"numCores":   4, "totalminMemoryMB":   30425, "ephemeralStorageGB":   72, "rank": 6}},
        "mem3_ssd1_x8":       {"internalName": "r3.2xlarge", "traits": {"numCores":   8, "totalminMemoryMB":   61187, "ephemeralStorageGB":  147, "rank": 9}},
        "mem3_ssd1_x16":      {"internalName": "r3.4xlarge", "traits": {"numCores":  16, "totalminMemoryMB":  122705, "ephemeralStorageGB":  297, "rank": 11}},
        "mem3_ssd1_x32":      {"internalName": "r3.8xlarge", "traits": {"numCores":  32, "totalminMemoryMB":  245751, "ephemeralStorageGB":  597, "rank": 13}}
}"""

  // Create an availble instance list based on a hard coded list
  private lazy val testDb: InstanceTypeDB = {
    val allInstances: Map[String, JsValue] = instanceList.parseJson.asJsObject.fields
    val db = allInstances.map {
      case (name, v) =>
        val fields: Map[String, JsValue] = v.asJsObject.fields
        val traits = fields("traits").asJsObject.fields
        val minMemoryMB = JsUtils.getInt(traits, "totalminMemoryMB")
        val diskGB = JsUtils.getInt(traits, "ephemeralStorageGB")
        val cpu = JsUtils.getInt(traits, "numCores")
        val priceRank = JsUtils.getInt(traits, "rank")
        val diskType = name.toLowerCase match {
          case s if s.contains("_ssd") => Some(DiskType.SSD)
          case s if s.contains("_hdd") => Some(DiskType.HDD)
          case _                       => None
        }
        name -> DxInstanceType(name,
                               minMemoryMB,
                               diskGB,
                               cpu,
                               gpu = false,
                               Vector(Constants.DefaultExecutionEnvironment),
                               diskType,
                               Some(priceRank))
    }
    InstanceTypeDB(db)
  }

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
    def makeString(s: String): TAT.Expr =
      TAT.ValueString(s, WdlTypes.T_String, Quoting.Double)(SourceLocation.empty)
    val rt = SeqMap(
        Runtime.DxInstanceTypeKey -> dxInstanceType.map(makeString),
        WdlRuntime.MemoryKey -> memory.map(makeString),
        WdlRuntime.DisksKey -> disks.map(makeString),
        WdlRuntime.CpuKey -> cpu.map(makeString),
        WdlRuntime.GpuKey -> gpu.map(b =>
          TAT.ValueBoolean(b, WdlTypes.T_Boolean)(SourceLocation.empty)
        )
    ).collect {
      case (key, Some(value)) => key -> value
    }
    Runtime(
        WdlVersion.V1,
        Some(TAT.RuntimeSection(rt)(SourceLocation.empty)),
        None,
        evaluator
    )
  }

  it should "Choose reasonable platform instance types" in {
    testDb.selectOptimal(InstanceTypeRequest.empty) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem1_ssd1_x2" =>
    }
    testDb.selectOptimal(
        InstanceTypeRequest(minMemoryMB = Some(3 * 1024), minDiskGB = Some(100), minCpu = Some(5))
    ) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem1_ssd1_x8" =>
    }
    testDb
      .selectOptimal(InstanceTypeRequest(minMemoryMB = Some(2 * 1024), minDiskGB = Some(20))) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem1_ssd1_x2" =>
    }
    testDb.selectOptimal(
        InstanceTypeRequest(minMemoryMB = Some(30 * 1024), minDiskGB = Some(128), minCpu = Some(8))
    ) should matchPattern {
      case Some(instanceType: DxInstanceType) if instanceType.name == "mem3_ssd1_x8" =>
    }

    // no instance with 1024 CPUs
    testDb.selectOptimal(InstanceTypeRequest(minCpu = Some(1024))) shouldBe None

    testDb
      .apply(
          createRuntime(None, Some("3 GB"), Some("local-disk 10 HDD"), Some("1"), None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x2")
    testDb
      .apply(
          createRuntime(None, Some("37 GB"), Some("local-disk 10 HDD"), Some("6"), None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x8")
    testDb
      .apply(
          createRuntime(None, Some("2 GB"), Some("local-disk 100 HDD"), None, None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x8")
    testDb
      .apply(
          createRuntime(None, Some("2.1GB"), Some("local-disk 100 HDD"), None, None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x8")

    testDb
      .apply(
          createRuntime(Some("mem3_ssd1_x8"), None, None, None, None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x8")

    testDb
      .apply(
          createRuntime(None, Some("235 GB"), Some("local-disk 550 HDD"), Some("32"), None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x32")
    testDb
      .apply(
          createRuntime(Some("mem3_ssd1_x32"), None, None, None, None).parseInstanceType
      )
      .name should equal("mem3_ssd1_x32")

    testDb
      .apply(
          createRuntime(None, None, None, Some("8"), None).parseInstanceType
      )
      .name should equal("mem1_ssd1_x8")
  }

  it should "catch parsing errors" in {
    // TODO: due to coercion, an int value is able to be coerced to a string, so this doesn't
    //  actually throw an exception at parse time, though it will throw an exception when the value
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
          minCpu = Some(1),
          os = Some(Constants.DefaultExecutionEnvironment)
      )

    createRuntime(None, Some("230MiB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(230),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          os = Some(Constants.DefaultExecutionEnvironment))

    createRuntime(None, Some("230GB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(
          minMemoryMB = Some(
              Math.ceil((230d * 1000d * 1000d * 1000d) / (1024d * 1024d)).toLong
          ),
          minDiskGB = Some(1),
          minCpu = Some(1),
          os = Some(Constants.DefaultExecutionEnvironment)
      )

    createRuntime(None, Some("230GiB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(230 * 1024),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          os = Some(Constants.DefaultExecutionEnvironment))

    createRuntime(None, Some("1000 TB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(
          minMemoryMB =
            Some(Math.ceil((1000d * 1000d * 1000d * 1000d * 1000d) / (1024d * 1024d)).toLong),
          minDiskGB = Some(1),
          minCpu = Some(1),
          os = Some(Constants.DefaultExecutionEnvironment)
      )

    createRuntime(None, Some("1000 TiB"), None, None, None).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(1000L * 1024L * 1024L),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          os = Some(Constants.DefaultExecutionEnvironment))

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
        minCpu = Some(1),
        os = Some(Constants.DefaultExecutionEnvironment)
    )
    createRuntime(None, None, None, Some("1.2"), None).parseInstanceType shouldBe InstanceTypeRequest(
        minMemoryMB = Some(2048),
        minDiskGB = Some(1),
        minCpu = Some(2),
        os = Some(Constants.DefaultExecutionEnvironment)
    )

    //    assertThrows[Exception] {
    //      createRuntime(None, None, None, Some(V_Boolean(false)), None)
    //    }

    // gpu
    createRuntime(None, Some("1000 TiB"), None, None, Some(true)).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(1000L * 1024L * 1024L),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          gpu = Some(true),
                          os = Some(Constants.DefaultExecutionEnvironment))

    createRuntime(None, None, None, None, Some(false)).parseInstanceType shouldBe
      InstanceTypeRequest(minMemoryMB = Some(2048),
                          minDiskGB = Some(1),
                          minCpu = Some(1),
                          gpu = Some(false),
                          os = Some(Constants.DefaultExecutionEnvironment))
  }

  it should "get required instance type" in {
    val db = InstanceType.createDb(Some(DxApi.get.project("project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq")))
    // query: memory=Some(2048) disk=Some(1) diskType=None cores=Some(5) gpu=None  os=None instancetype=None
    val runtime = Runtime(
        WdlVersion.V1,
        Some(
            TAT.RuntimeSection(
                SeqMap("cpu" -> TAT.ValueInt(5, WdlTypes.T_Int)(SourceLocation.empty))
            )(
                SourceLocation.empty
            )
        ),
        None,
        evaluator,
        None,
        None
    )
    val request = runtime.parseInstanceType
    // TODO: we have requested that v2 instance types be enabled for the staging org - this test
    //  will break when that happens as the selected instance type should be mem1_ssd1_v2_x8
    db.apply(request).name shouldBe "mem1_ssd1_x8"
  }
}
