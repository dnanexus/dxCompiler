package dx.executor

import dx.Assumptions.{isLoggedIn, toolkitCallable}
import dx.api.{DxApi, DxFileDescribe, Field}
import dx.util.{FileUtils, Logger}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Files
import scala.util.Random

class FileUploaderTest extends AnyFlatSpec with Matchers with BeforeAndAfterAll {
  assume(isLoggedIn)
  assume(toolkitCallable)
  private val logger = Logger.Quiet
  private val dxApi = DxApi()(logger)
  val testProject = "dxCompiler_playground"

  private val dxTestProject =
    try {
      dxApi.resolveProject(testProject)
    } catch {
      case _: Exception =>
        throw new Exception(
            s"""|Could not find project ${testProject}, you probably need to be logged into
                |the platform""".stripMargin
        )
    }

  private val username = dxApi.whoami()
  private val uploadPath = s"unit_tests/${username}/test_upload"
  private val testDir = Files.createTempDirectory("test")
  private val random = new Random(42)
  private val files = Iterator
    .range(0, 5)
    .map { i =>
      val path = testDir.resolve(s"file_${i}.txt")
      val length = random.nextInt(1024 * 1024)
      val content = random.nextString(length)
      FileUtils.writeFileContent(path, content)
      (path, length)
    }
    .toMap
  private val fileTags = Set("A", "B")
  private val fileProperties = Map("name" -> "Joe", "age" -> "42")

  override protected def afterAll(): Unit = {
    dxTestProject.removeFolder(s"/${uploadPath}/", recurse = true, force = true)
  }

  it should "upload files in serial" in {
    val uploader = SerialFileUploader(dxApi)
    val dest = s"${dxTestProject.id}:/${uploadPath}/serial/"
    val uploads = files.keys.map { path =>
      FileUpload(path, Some(dest), fileTags, fileProperties)
    }.toSet
    val results = uploader.upload(uploads, wait = true)
    results.size shouldBe 5
    results.foreach {
      case (path, dxFile) =>
        val desc: DxFileDescribe = dxFile.describe(Set(Field.Tags, Field.Properties))
        files(path) shouldBe desc.size
        desc.tags shouldBe Some(fileTags)
        desc.properties shouldBe Some(fileProperties)
    }
  }

  it should "upload files in parallel" in {
    val uploader = ParallelFileUploader(maxConcurrent = 3, dxApi = dxApi)
    val dest = s"${dxTestProject.id}:/${uploadPath}/parallel/"
    val uploads = files.keys.map { path =>
      FileUpload(path, Some(dest), fileTags, fileProperties)
    }.toSet
    val results = uploader.upload(uploads, wait = true)
    results.size shouldBe 5
    results.foreach {
      case (path, dxFile) =>
        val desc: DxFileDescribe = dxFile.describe(Set(Field.Tags, Field.Properties))
        files(path) shouldBe desc.size
        desc.tags shouldBe Some(fileTags)
        desc.properties shouldBe Some(fileProperties)
    }
  }
}
