package dx.core.io

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import java.nio.charset.StandardCharsets

/**
  * Unit tests for AuthenticatedHttpFileAccessProtocol.
  * These tests verify the protocol's behavior without making actual HTTP requests.
  */
class AuthenticatedHttpFileAccessProtocolTest extends AnyFlatSpec with Matchers {

  private val testToken = "test-token-12345"
  private val otherToken = "other-token-67890"

  "AuthenticatedHttpFileAccessProtocol" should "authenticate to configured domains" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      domainTokens = Map("github.com" -> testToken, "raw.githubusercontent.com" -> testToken)
    )

    val source = protocol.resolve("https://raw.githubusercontent.com/org/repo/main/file.wdl")
    source.token shouldBe Some(testToken)
  }

  it should "not authenticate to unconfigured domains" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      domainTokens = Map("github.com" -> testToken)
    )

    val source = protocol.resolve("https://example.com/file.wdl")
    source.token shouldBe None
  }

  it should "use different tokens for different domains" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      domainTokens = Map(
        "raw.githubusercontent.com" -> testToken,
        "gitlab.com" -> otherToken
      )
    )

    val ghSource = protocol.resolve("https://raw.githubusercontent.com/org/repo/main/file.wdl")
    ghSource.token shouldBe Some(testToken)

    val glSource = protocol.resolve("https://gitlab.com/org/repo/file.wdl")
    glSource.token shouldBe Some(otherToken)
  }

  it should "work without any tokens (backward compatible)" in {
    val protocol = AuthenticatedHttpFileAccessProtocol()

    val source = protocol.resolve("https://github.com/org/repo/file.wdl")
    source.token shouldBe None
  }

  it should "be case-insensitive for domain matching" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      domainTokens = Map("github.com" -> testToken)
    )

    val source = protocol.resolve("https://GitHub.COM/org/repo/file.wdl")
    source.token shouldBe Some(testToken)
  }

  it should "support directory resolution with auth" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      domainTokens = Map("github.com" -> testToken)
    )

    val dirSource = protocol.resolveDirectory("https://github.com/org/repo/archive.tar.gz")
    dirSource.isDirectory shouldBe true
    dirSource.token shouldBe Some(testToken)
  }

  it should "handle HTTP and HTTPS schemes" in {
    val protocol = AuthenticatedHttpFileAccessProtocol()
    protocol.schemes should contain("http")
    protocol.schemes should contain("https")
  }

  it should "support directories" in {
    val protocol = AuthenticatedHttpFileAccessProtocol()
    protocol.supportsDirectories shouldBe true
  }

  it should "handle empty domain tokens map" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(domainTokens = Map.empty)

    val source = protocol.resolve("https://github.com/org/repo/file.wdl")
    source.token shouldBe None
  }

  "parseTokens" should "parse semicolon-separated domain:token pairs" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "raw.githubusercontent.com:ghp_abc123;gitlab.com:glpat-xyz789"
    )
    tokens shouldBe Map(
      "raw.githubusercontent.com" -> "ghp_abc123",
      "gitlab.com" -> "glpat-xyz789"
    )
  }

  it should "handle a single domain:token pair" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "raw.githubusercontent.com:ghp_abc123"
    )
    tokens shouldBe Map("raw.githubusercontent.com" -> "ghp_abc123")
  }

  it should "handle extra whitespace" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "  github.com : token1 ;  gitlab.com : token2  "
    )
    tokens shouldBe Map("github.com" -> "token1", "gitlab.com" -> "token2")
  }

  it should "handle trailing semicolons" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "github.com:token1;"
    )
    tokens shouldBe Map("github.com" -> "token1")
  }

  it should "skip malformed entries without colons" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "github.com:token1;badentry;gitlab.com:token2"
    )
    tokens shouldBe Map("github.com" -> "token1", "gitlab.com" -> "token2")
  }

  it should "split only on first colon (tokens may contain colons)" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "github.com:token:with:colons"
    )
    tokens shouldBe Map("github.com" -> "token:with:colons")
  }

  it should "lowercase domain names" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      "GitHub.COM:mytoken"
    )
    tokens shouldBe Map("github.com" -> "mytoken")
  }

  it should "return empty map for empty string" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens("")
    tokens shouldBe Map.empty
  }

  it should "skip entries with empty domain or token" in {
    val tokens = AuthenticatedHttpFileAccessProtocol.parseTokens(
      ":token1;github.com:"
    )
    tokens shouldBe Map.empty
  }

  "AuthenticatedHttpFileSource" should "resolve relative paths with token" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/"),
      StandardCharsets.UTF_8,
      isDirectory = true,
      token = Some(testToken)
    )("https://github.com/org/repo/main/")

    val resolved = source.resolve("subdir/file.wdl")
    resolved.address should include("subdir/file.wdl")
    resolved.token shouldBe Some(testToken)
  }

  it should "propagate token to resolved files" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/"),
      StandardCharsets.UTF_8,
      isDirectory = true,
      token = Some(testToken)
    )("https://github.com/org/repo/main/")

    val resolved = source.resolve("another_file.wdl")
    resolved.token shouldBe Some(testToken)
  }

  it should "propagate token to resolved directories" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/"),
      StandardCharsets.UTF_8,
      isDirectory = true,
      token = Some(testToken)
    )("https://github.com/org/repo/main/")

    val resolved = source.resolveDirectory("subdir")
    resolved.isDirectory shouldBe true
    resolved.token shouldBe Some(testToken)
  }

  it should "get parent directory with token" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/file.wdl"),
      StandardCharsets.UTF_8,
      isDirectory = false,
      token = Some(testToken)
    )("https://github.com/org/repo/main/file.wdl")

    val parent = source.getParent
    parent shouldBe defined
    parent.get.isDirectory shouldBe true
    parent.get.token shouldBe Some(testToken)
  }

  it should "extract name from URI" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/file.wdl"),
      StandardCharsets.UTF_8,
      isDirectory = false,
      token = None
    )("https://github.com/org/repo/main/file.wdl")

    source.name shouldBe "file.wdl"
  }

  it should "extract folder from URI" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/file.wdl"),
      StandardCharsets.UTF_8,
      isDirectory = false,
      token = None
    )("https://github.com/org/repo/main/file.wdl")

    source.folder shouldBe "/org/repo/main"
  }

  it should "not be listable" in {
    val source = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/"),
      StandardCharsets.UTF_8,
      isDirectory = true,
      token = None
    )("https://github.com/org/repo/main/")

    source.isListable shouldBe false
  }

  it should "relativize paths correctly" in {
    val dirSource = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/"),
      StandardCharsets.UTF_8,
      isDirectory = true,
      token = None
    )("https://github.com/org/repo/main/")

    val fileSource = AuthenticatedHttpFileSource(
      java.net.URI.create("https://github.com/org/repo/main/subdir/file.wdl"),
      StandardCharsets.UTF_8,
      isDirectory = false,
      token = None
    )("https://github.com/org/repo/main/subdir/file.wdl")

    val relativePath = dirSource.relativize(fileSource)
    relativePath shouldBe "subdir/file.wdl"
  }
}
