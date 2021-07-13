package dx.executor

import dx.Assumptions.{isLoggedIn, toolkitCallable}
import dx.api.DxApi
import dx.util.Logger
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class FileUploaderTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)
  assume(toolkitCallable)
  private val logger = Logger.Quiet
  private val dxApi = DxApi()(logger)
  val testProject = "dxCompiler_playground"

  private lazy val dxTestProject =
    try {
      dxApi.resolveProject(testProject)
    } catch {
      case _: Exception =>
        throw new Exception(
            s"""|Could not find project ${testProject}, you probably need to be logged into
                |the platform""".stripMargin
        )
    }

  private lazy val username = dxApi.whoami()
  private lazy val unitTestsPath = s"unit_tests/${username}"
  it should "upload files in serial" in {}

  it should "upload files in parallel" in {}
}
