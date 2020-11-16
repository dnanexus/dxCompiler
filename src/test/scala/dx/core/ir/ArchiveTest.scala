package dx.core.ir

import java.nio.file.Files

import dx.util.FileUtils
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class ArchiveTest extends AnyFlatSpec with Matchers {
  it should "create a squasfs" in {
    val tmpDir = Files.createTempDirectory("archive")
    val subDir = tmpDir.resolve("sub")
    val file = subDir.resolve("file.txt")
    FileUtils.writeFileContent(file, "some text")
    val archiveFile = tmpDir.resolve("test.img")
    val fs = SquashFs(archiveFile)()
    fs.isMounted shouldBe false
    fs.append(file, Some(tmpDir), removeOriginal = true)
    Files.exists(file) shouldBe false
    fs.mount()
    fs.isMounted shouldBe true
    val fileInFs = fs.resolve("sub/file.txt")
    FileUtils.copyFile(fileInFs, file)
    Files.exists(file) shouldBe true
    FileUtils.readFileContent(file) shouldBe "some text"
  }

  it should "create an archive" in {}
}
