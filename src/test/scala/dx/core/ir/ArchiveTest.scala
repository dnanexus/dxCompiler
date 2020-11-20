package dx.core.ir

import java.nio.file.{Files, Paths}

import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.FileUtils
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class ArchiveTest extends AnyFlatSpec with Matchers {
  it should "create a squasfs" in {
    val tmpDir = Files.createTempDirectory("archive")
    val file = tmpDir.resolve("file.txt")
    FileUtils.writeFileContent(file, "some text")
    val archiveFile = tmpDir.resolve("test.img")
    val fs = SquashFs(archiveFile)()
    fs.isMounted shouldBe false
    fs.append(file, Some(Paths.get("sub")), removeSource = true)
    Files.exists(file) shouldBe false
    fs.mount()
    fs.isMounted shouldBe true
    try {
      val fileInFs = fs.resolve("sub/file.txt")
      FileUtils.copyFile(fileInFs, file)
      Files.exists(file) shouldBe true
      FileUtils.readFileContent(file) shouldBe "some text"
    } finally {
      fs.unmount()
    }
    fs.isMounted shouldBe false
  }

  it should "create an archive" in {
    val tmpDir = Files.createTempDirectory("archive")
    val subDir1 = tmpDir.resolve("sub1")
    val file1 = subDir1.resolve("file1.txt")
    val file2 = subDir1.resolve("file2.txt")
    val subDir2 = tmpDir.resolve("sub2")
    val file3 = subDir2.resolve("file3.txt")
    FileUtils.writeFileContent(file1, "file1")
    FileUtils.writeFileContent(file2, "file2")
    FileUtils.writeFileContent(file3, "file3")
    val structType = TSchema("Files", Map("files" -> TArray(TFile)))
    val value = VHash(
        Map(
            "files" -> VArray(
                Vector(VFile(file1.toString), VFile(file2.toString), VFile(file3.toString))
            )
        )
    )
    val unpacked = LocalizedArchive(structType, value)(parentDir = Some(tmpDir), name = Some("foo"))
    unpacked.localized shouldBe true
    unpacked.isOpen shouldBe false
    val packed = unpacked.pack()
    packed.isOpen shouldBe false
    packed.localized shouldBe false
    packed.irType shouldBe structType
    packed.irValue shouldBe VHash(
        Map(
            "files" -> VArray(
                Vector(
                    VFile(Paths.get("sub1").resolve("file1.txt").toString),
                    VFile(Paths.get("sub1").resolve("file2.txt").toString),
                    VFile(Paths.get("sub2").resolve("file3.txt").toString)
                )
            )
        )
    )
    val packed2 = PackedArchive(packed.path)(typeAliases = Some(Map.empty))
    packed2.isOpen shouldBe false
    packed2.localized shouldBe false
    packed2.irType shouldBe structType
    packed2.irValue shouldBe VHash(
        Map(
            "files" -> VArray(
                Vector(
                    VFile(Paths.get("sub1").resolve("file1.txt").toString),
                    VFile(Paths.get("sub1").resolve("file2.txt").toString),
                    VFile(Paths.get("sub2").resolve("file3.txt").toString)
                )
            )
        )
    )
    val (localized, _) = packed2.localize()
    try {
      localized.isOpen shouldBe true
      localized.irType shouldBe structType
      val localizedFiles = localized.irValue match {
        case VHash(m) =>
          m.get("files") match {
            case Some(VArray(files)) =>
              files.map {
                case VFile(path) => path
                case _           => throw new Exception("invalid localized value")
              }
            case _ => throw new Exception("invalid localized value")
          }
        case _ => throw new Exception("invalid localized value")
      }
      val localizedFile1 = Paths.get(localizedFiles(0))
      localizedFile1.getFileName.toString shouldBe "file1.txt"
      FileUtils.readFileContent(localizedFile1) shouldBe "file1"
      val localizedFile2 = Paths.get(localizedFiles(1))
      localizedFile2.getFileName.toString shouldBe "file2.txt"
      FileUtils.readFileContent(localizedFile2) shouldBe "file2"
      val localizedFile3 = Paths.get(localizedFiles(2))
      localizedFile3.getFileName.toString shouldBe "file3.txt"
      FileUtils.readFileContent(localizedFile3) shouldBe "file3"
    } finally {
      localized.close()
    }
  }
}
