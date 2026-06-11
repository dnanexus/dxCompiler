package dx.core.languages.wdl

import dx.util.{FileUtils, Logger}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class WdlImportHttpAuthTest extends AnyFlatSpec with Matchers {

  it should "expose the dxCompiler-prefixed env variable name" in {
    WdlImportHttpAuth.TokensEnvVar shouldBe "DXCOMPILER_WDL_IMPORT_BEARER_TOKENS"
  }

  it should "return a protocol with no tokens when DXCOMPILER_WDL_IMPORT_BEARER_TOKENS is unset" in {
    // We can't portably mutate process env vars on JVM 11, so we only verify
    // the unset path. The set path is logically equivalent to parseAuthTokens(value)
    // followed by direct construction, both covered by AuthenticatedHttpFileAccessProtocolTest.
    assume(sys.env.get(WdlImportHttpAuth.TokensEnvVar).isEmpty,
           s"${WdlImportHttpAuth.TokensEnvVar} is set in the test environment; skipping unset-path test")
    val protocol = WdlImportHttpAuth.fromEnvironment()
    protocol.domainBearerTokens shouldBe empty
    protocol.resolve("https://raw.githubusercontent.com/x.wdl").credentials shouldBe None
  }

  it should "configure the protocol with the unauthorized hint" in {
    assume(sys.env.get(WdlImportHttpAuth.TokensEnvVar).isEmpty)
    val protocol = WdlImportHttpAuth.fromEnvironment()
    protocol.unauthorizedHint shouldBe Some(WdlImportHttpAuth.UnauthorizedHint)
  }

  it should "forward the unauthorized hint to constructed file sources so the 401 message can guide users" in {
    assume(sys.env.get(WdlImportHttpAuth.TokensEnvVar).isEmpty)
    val protocol = WdlImportHttpAuth.fromEnvironment()
    protocol
      .resolve("https://raw.githubusercontent.com/x.wdl")
      .unauthorizedHint shouldBe Some(WdlImportHttpAuth.UnauthorizedHint)
    protocol
      .resolveDirectory("https://raw.githubusercontent.com/dir/")
      .unauthorizedHint shouldBe Some(WdlImportHttpAuth.UnauthorizedHint)
  }

  it should "use the default encoding and a Quiet logger when called with no args" in {
    assume(sys.env.get(WdlImportHttpAuth.TokensEnvVar).isEmpty)
    val protocol = WdlImportHttpAuth.fromEnvironment()
    protocol.encoding shouldBe FileUtils.DefaultEncoding
    protocol.logger shouldBe Logger.Quiet
  }
}
