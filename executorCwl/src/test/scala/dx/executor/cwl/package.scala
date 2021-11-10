package dx.executor.cwl

import dx.api.DxApi
import dx.util.SysUtils
import org.scalatest.Tag

// test that requires being logged into DNAnexus
class DxTag(name: String = "dx") extends Tag(name)

object Tags {
  // a test that calls a DNAnexus API method (and thus requires being
  // logged in), but doesn't create any objects
  object ApiTest extends DxTag("dxApi")
  // test that builds native applets/workflows
  object NativeTest extends DxTag("native")
  // test that rqeuires being logged into a DNAnexus production account
  object ProdTest extends DxTag("prod")
  // marker for an edge case
  object EdgeTest extends Tag("edge")
}

object Assumptions {
  private val dxApi = DxApi.get

  /**
    * Tests that the user is logged in.
    * @return
    */
  def isLoggedIn: Boolean = dxApi.isLoggedIn

  /**
    * Tests that the dx toolkit is installed and in the path.
    * @return
    */
  lazy val toolkitCallable: Boolean = {
    try {
      val (retcode, _, _) = SysUtils.execCommand("dx whoami")
      retcode == 0
    } catch {
      case _: Throwable => false
    }
  }

  /**
    * Tests that dxda is installed and in the path.
    */
  lazy val dxdaCallable: Boolean = {
    try {
      SysUtils.execCommand("dx-download-agent version")
      true
    } catch {
      case _: Throwable => false
    }
  }

  lazy val cwltoolCallable: Boolean = {
    val (retcode, stdout, stderr) =
      SysUtils.execCommand("cwltool --version", exceptionOnFailure = false)
    // TODO: pull the min cwl version from the bundled_dependencies.json file
    if (retcode != 0) {
      throw new Exception(
          f"cwltool --version failed with return code ${retcode}\nstdout:\n${stdout}\nstderr:\n${stderr}"
      )
    }
    val installedVersion = stdout.trim.split('.').last.toLong
    if (installedVersion < 20210628163208L) {
      throw new Exception(f"wrong cwltool version installed ${installedVersion}")
    }
    true
  }
}
