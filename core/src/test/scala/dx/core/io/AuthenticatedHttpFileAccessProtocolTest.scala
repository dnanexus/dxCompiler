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

  "AuthenticatedHttpFileAccessProtocol" should "read token from constructor" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set("github.com")
    )
    protocol.token shouldBe Some(testToken)
  }

  it should "only authenticate to allowed domains" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set("github.com", "raw.githubusercontent.com")
    )

    // GitHub should get auth
    val githubSource = protocol.resolve("https://raw.githubusercontent.com/org/repo/main/file.wdl")
    githubSource.token shouldBe Some(testToken)

    // Other domains should NOT get auth
    val otherSource = protocol.resolve("https://example.com/file.wdl")
    otherSource.token shouldBe None
  }

  it should "work without a token (backward compatible)" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = None,
      allowedDomains = Set("github.com")
    )

    val source = protocol.resolve("https://github.com/org/repo/file.wdl")
    source.token shouldBe None
  }

  it should "be case-insensitive for domain matching" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set("github.com")
    )

    val source = protocol.resolve("https://GitHub.COM/org/repo/file.wdl")
    source.token shouldBe Some(testToken)
  }

  it should "support directory resolution" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set("github.com")
    )

    val dirSource = protocol.resolveDirectory("https://github.com/org/repo/archive.tar.gz")
    dirSource.isDirectory shouldBe true
    dirSource.token shouldBe Some(testToken)
  }

  it should "handle HTTP scheme" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set("example.com")
    )

    protocol.schemes should contain("http")
    protocol.schemes should contain("https")
  }

  it should "support directories" in {
    val protocol = AuthenticatedHttpFileAccessProtocol()
    protocol.supportsDirectories shouldBe true
  }

  it should "not send token to unlisted domains even with token configured" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set("github.com")
    )

    val source = protocol.resolve("https://gitlab.com/org/repo/file.wdl")
    source.token shouldBe None
  }

  it should "handle empty allowed domains set" in {
    val protocol = AuthenticatedHttpFileAccessProtocol(
      token = Some(testToken),
      allowedDomains = Set.empty
    )

    val source = protocol.resolve("https://github.com/org/repo/file.wdl")
    source.token shouldBe None
  }

  "AuthenticatedHttpFileAccessProtocol.defaultDomains" should "include github.com" in {
    AuthenticatedHttpFileAccessProtocol.defaultDomains should contain("github.com")
  }

  it should "include raw.githubusercontent.com" in {
    AuthenticatedHttpFileAccessProtocol.defaultDomains should contain("raw.githubusercontent.com")
  }

  "Domain parsing" should "parse comma-separated domains correctly" in {
    val domainsString = "gitlab.com, bitbucket.org, custom.example.com"
    val parsed = domainsString.split(",").map(_.trim.toLowerCase).filter(_.nonEmpty).toSet

    parsed should contain("gitlab.com")
    parsed should contain("bitbucket.org")
    parsed should contain("custom.example.com")
    parsed.size shouldBe 3
  }

  it should "handle extra whitespace in domain list" in {
    val domainsString = "  github.com  ,  gitlab.com  ,  "
    val parsed = domainsString.split(",").map(_.trim.toLowerCase).filter(_.nonEmpty).toSet

    parsed should contain("github.com")
    parsed should contain("gitlab.com")
    parsed.size shouldBe 2
  }

  "AuthenticatedHttpFileSource" should "resolve relative paths" in {
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
