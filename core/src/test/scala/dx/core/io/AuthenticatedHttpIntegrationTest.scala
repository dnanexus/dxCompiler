package dx.core.io

import com.sun.net.httpserver.{HttpExchange, HttpHandler, HttpServer}
import dx.util.{FileSourceResolver, LocalFileAccessProtocol, Logger}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.net.InetSocketAddress
import java.nio.charset.StandardCharsets

/**
  * Integration tests using an embedded HTTP server to verify
  * authenticated imports work correctly end-to-end.
  */
class AuthenticatedHttpIntegrationTest extends AnyFlatSpec with Matchers with BeforeAndAfterAll {

  private var server: HttpServer = _
  private var serverPort: Int = _
  private val testToken = "test-bearer-token-xyz"

  private val wdlContent =
    """version 1.0
      |
      |task hello {
      |  command { echo "Hello, World!" }
      |  output { String message = read_string(stdout()) }
      |}
      |""".stripMargin

  override def beforeAll(): Unit = {
    super.beforeAll()

    // Create embedded HTTP server
    server = HttpServer.create(new InetSocketAddress(0), 0)
    serverPort = server.getAddress.getPort

    // Public endpoint - no auth required
    server.createContext("/public/file.wdl", new HttpHandler {
      override def handle(exchange: HttpExchange): Unit = {
        val response = wdlContent.getBytes(StandardCharsets.UTF_8)
        exchange.sendResponseHeaders(200, response.length)
        val os = exchange.getResponseBody
        os.write(response)
        os.close()
      }
    })

    // Private endpoint - requires Bearer token
    server.createContext("/private/file.wdl", new HttpHandler {
      override def handle(exchange: HttpExchange): Unit = {
        val authHeader = exchange.getRequestHeaders.getFirst("Authorization")

        if (authHeader == s"Bearer $testToken") {
          val response = wdlContent.getBytes(StandardCharsets.UTF_8)
          exchange.sendResponseHeaders(200, response.length)
          val os = exchange.getResponseBody
          os.write(response)
          os.close()
        } else if (authHeader == null) {
          exchange.sendResponseHeaders(401, -1)
          exchange.close()
        } else {
          exchange.sendResponseHeaders(403, -1)
          exchange.close()
        }
      }
    })

    // Endpoint that checks for token and returns different content
    server.createContext("/conditional/file.wdl", new HttpHandler {
      override def handle(exchange: HttpExchange): Unit = {
        val authHeader = exchange.getRequestHeaders.getFirst("Authorization")
        val content = if (authHeader == s"Bearer $testToken") {
          "authenticated content"
        } else {
          "public content"
        }
        val response = content.getBytes(StandardCharsets.UTF_8)
        exchange.sendResponseHeaders(200, response.length)
        val os = exchange.getResponseBody
        os.write(response)
        os.close()
      }
    })

    // HEAD endpoint for exists check
    server.createContext("/head-test/file.wdl", new HttpHandler {
      override def handle(exchange: HttpExchange): Unit = {
        val authHeader = exchange.getRequestHeaders.getFirst("Authorization")

        if (authHeader == s"Bearer $testToken") {
          exchange.sendResponseHeaders(200, -1)
        } else if (authHeader == null) {
          exchange.sendResponseHeaders(401, -1)
        } else {
          exchange.sendResponseHeaders(403, -1)
        }
        exchange.close()
      }
    })

    server.setExecutor(null)
    server.start()
  }

  override def afterAll(): Unit = {
    if (server != null) {
      server.stop(0)
    }
    super.afterAll()
  }

  private def createProtocol(token: Option[String], allowedDomains: Set[String]): AuthenticatedHttpFileAccessProtocol = {
    AuthenticatedHttpFileAccessProtocol(
      token = token,
      allowedDomains = allowedDomains,
      logger = Logger.Quiet
    )
  }

  "Authenticated HTTP imports" should "access public endpoints without token" in {
    val protocol = createProtocol(None, Set.empty)
    val source = protocol.resolve(s"http://localhost:$serverPort/public/file.wdl")

    source.readString should include("version 1.0")
  }

  it should "access private endpoints with valid token" in {
    val protocol = createProtocol(Some(testToken), Set("localhost"))
    val source = protocol.resolve(s"http://localhost:$serverPort/private/file.wdl")

    source.readString should include("version 1.0")
  }

  it should "fail on private endpoints without token" in {
    val protocol = createProtocol(None, Set.empty)
    val source = protocol.resolve(s"http://localhost:$serverPort/private/file.wdl")

    val exception = intercept[Exception] {
      source.readString
    }
    exception.getMessage should include("401")
    exception.getMessage should include("WDL_IMPORT_TOKEN")
  }

  it should "fail on private endpoints with wrong token" in {
    val protocol = createProtocol(Some("wrong-token"), Set("localhost"))
    val source = protocol.resolve(s"http://localhost:$serverPort/private/file.wdl")

    val exception = intercept[Exception] {
      source.readString
    }
    exception.getMessage should include("403")
  }

  it should "not send token to non-allowed domains" in {
    // Token configured but localhost not in allowed domains
    val protocol = createProtocol(Some(testToken), Set("github.com"))
    val source = protocol.resolve(s"http://localhost:$serverPort/conditional/file.wdl")

    // Should get public content since token wasn't sent
    source.readString shouldBe "public content"
  }

  it should "send token only to allowed domains" in {
    val protocol = createProtocol(Some(testToken), Set("localhost"))
    val source = protocol.resolve(s"http://localhost:$serverPort/conditional/file.wdl")

    // Should get authenticated content since token was sent
    source.readString shouldBe "authenticated content"
  }

  it should "read file content correctly" in {
    val protocol = createProtocol(Some(testToken), Set("localhost"))
    val source = protocol.resolve(s"http://localhost:$serverPort/private/file.wdl")

    val content = source.readString
    content should include("version 1.0")
    content should include("task hello")
    content should include("Hello, World!")
  }

  it should "check exists with authentication" in {
    val protocol = createProtocol(Some(testToken), Set("localhost"))
    val source = protocol.resolve(s"http://localhost:$serverPort/head-test/file.wdl")

    source.exists shouldBe true
  }

  it should "fail exists check without required token" in {
    val protocol = createProtocol(None, Set.empty)
    val source = protocol.resolve(s"http://localhost:$serverPort/head-test/file.wdl")

    val exception = intercept[Exception] {
      source.exists
    }
    exception.getMessage should include("401")
  }

  it should "fail exists check with wrong token" in {
    val protocol = createProtocol(Some("wrong-token"), Set("localhost"))
    val source = protocol.resolve(s"http://localhost:$serverPort/head-test/file.wdl")

    val exception = intercept[Exception] {
      source.exists
    }
    exception.getMessage should include("403")
  }

  it should "work with FileSourceResolver" in {
    val httpProtocol = createProtocol(Some(testToken), Set("localhost"))
    val resolver = FileSourceResolver(Vector(
      LocalFileAccessProtocol(),
      httpProtocol
    ))

    val source = resolver.resolve(s"http://localhost:$serverPort/private/file.wdl")
    source.readString should include("version 1.0")
  }

  it should "resolve relative imports with authentication" in {
    val protocol = createProtocol(Some(testToken), Set("localhost"))
    val baseSource = protocol.resolveDirectory(s"http://localhost:$serverPort/private/")

    val resolvedSource = baseSource.resolve("file.wdl")
    resolvedSource.token shouldBe Some(testToken)
  }
}
