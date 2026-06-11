package dx.core.languages.wdl

import com.sun.net.httpserver.{HttpExchange, HttpHandler, HttpServer}
import dx.api.DxApi
import dx.util.{
  AuthenticatedHttpFileAccessProtocol,
  FileSourceResolver,
  LocalFileAccessProtocol,
  Logger,
  StringFileNode
}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.net.InetSocketAddress
import java.nio.charset.StandardCharsets

/**
  * End-to-end test for authenticated WDL imports.
  *
  * Spins up a local HTTP server that guards a tiny WDL document behind a
  * static Bearer token, then exercises three configurations:
  *   - no token configured  -> 401
  *   - wrong token          -> 403
  *   - correct token        -> success (bytes round-trip and a main.wdl that
  *                                       imports the protected doc parses).
  */
class WdlImportHttpAuthIntegrationTest extends AnyFlatSpec with Matchers with BeforeAndAfterAll {

  private val ExpectedToken = "secret-token-abc"

  private val ImportedWdl: String =
    """version 1.0
      |
      |task hello {
      |  input {
      |    String name
      |  }
      |  command <<<
      |    echo "Hello, ~{name}"
      |  >>>
      |  output {
      |    String greeting = read_string(stdout())
      |  }
      |}
      |""".stripMargin

  private def importedWdlBytes: Array[Byte] = ImportedWdl.getBytes(StandardCharsets.UTF_8)

  private val protectedHandler: HttpHandler = (exchange: HttpExchange) => {
    val auth = Option(exchange.getRequestHeaders.getFirst("Authorization"))
    if (auth.contains(s"Bearer $ExpectedToken")) {
      val body = importedWdlBytes
      if (exchange.getRequestMethod == "HEAD") {
        exchange.sendResponseHeaders(200, -1)
        exchange.close()
      } else {
        exchange.sendResponseHeaders(200, body.length.toLong)
        val os = exchange.getResponseBody
        try os.write(body)
        finally os.close()
      }
    } else {
      // 401 when no credentials are supplied, 403 when credentials are supplied
      // but do not match.
      val status = if (auth.isEmpty) 401 else 403
      exchange.sendResponseHeaders(status, -1)
      exchange.close()
    }
  }

  private var server: HttpServer = _
  private var serverHost: String = _
  private var serverPort: Int = _

  override def beforeAll(): Unit = {
    server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0)
    server.createContext("/imported.wdl", protectedHandler)
    server.setExecutor(null)
    server.start()
    serverHost = server.getAddress.getHostString
    serverPort = server.getAddress.getPort
  }

  override def afterAll(): Unit = {
    if (server != null) server.stop(0)
  }

  private def importedUrl: String = s"http://${serverHost}:${serverPort}/imported.wdl"

  /** Build a resolver whose http/https handler is the authenticated protocol
    * with the given token map. The local protocol is preserved so any local
    * fallback still works. */
  private def resolverWith(tokens: Map[String, String]): FileSourceResolver = {
    val authProtocol = AuthenticatedHttpFileAccessProtocol(
        domainBearerTokens = tokens,
        unauthorizedHint = Some(WdlImportHttpAuth.UnauthorizedHint)
    )
    FileSourceResolver(
        Vector(
            LocalFileAccessProtocol(),
            authProtocol
        )
    )
  }

  // --- Direct resolver layer: fetch the protected WDL bytes ------------

  it should "fail with HTTP 401 and surface the env var hint when no token is configured" in {
    val resolver = resolverWith(Map.empty)
    val fs = resolver.resolve(importedUrl)
    val thrown = the[Exception] thrownBy fs.readBytes
    thrown.getMessage should include("HTTP 401 Unauthorized")
    thrown.getMessage should include(WdlImportHttpAuth.TokensEnvVar)
  }

  it should "fail with HTTP 403 when the configured token is wrong" in {
    val resolver = resolverWith(Map(serverHost -> "not-the-right-token"))
    val fs = resolver.resolve(importedUrl)
    val thrown = the[Exception] thrownBy fs.readBytes
    thrown.getMessage should include("HTTP 403 Forbidden")
  }

  it should "fetch the protected WDL bytes when the configured token is correct" in {
    val resolver = resolverWith(Map(serverHost -> ExpectedToken))
    val fs = resolver.resolve(importedUrl)
    new String(fs.readBytes, StandardCharsets.UTF_8) shouldBe ImportedWdl
  }

  it should "only attach the bearer token to hosts present in the token map" in {
    // Wrong-host entry must not cause the right host to be authenticated.
    val resolver = resolverWith(Map("some.other.host" -> ExpectedToken))
    val fs = resolver.resolve(importedUrl)
    val thrown = the[Exception] thrownBy fs.readBytes
    thrown.getMessage should include("HTTP 401 Unauthorized")
  }

  // --- WDL parser layer: a workflow that imports the protected document ---

  private def mainWdlImporting(remoteUrl: String): String =
    s"""version 1.0
       |
       |import "$remoteUrl" as lib
       |
       |workflow main {
       |  input {
       |    String name
       |  }
       |  call lib.hello { input: name = name }
       |  output {
       |    String greeting = hello.greeting
       |  }
       |}
       |""".stripMargin

  private def parseMain(resolver: FileSourceResolver) = {
    VersionSupport.fromSource(
        StringFileNode(mainWdlImporting(importedUrl)),
        WdlOptions.default,
        resolver,
        DxApi()(Logger.Quiet),
        Logger.Quiet
    )
  }

  it should "fail to parse a main.wdl that imports the protected doc when no token is set" in {
    val thrown = the[Exception] thrownBy parseMain(resolverWith(Map.empty))
    // The 401 surfaces somewhere in the cause chain; assert via the full chain string.
    val messages = Iterator
      .iterate(Option(thrown: Throwable))(_.flatMap(t => Option(t.getCause)))
      .takeWhile(_.isDefined)
      .flatMap(_.map(_.getMessage))
      .mkString(" | ")
    messages should (include("401") or include("Unauthorized"))
  }

  it should "fail to parse a main.wdl that imports the protected doc when the token is wrong" in {
    val thrown = the[Exception] thrownBy parseMain(
        resolverWith(Map(serverHost -> "wrong"))
    )
    val messages = Iterator
      .iterate(Option(thrown: Throwable))(_.flatMap(t => Option(t.getCause)))
      .takeWhile(_.isDefined)
      .flatMap(_.map(_.getMessage))
      .mkString(" | ")
    messages should (include("401") or include("Unauthorized"))
  }

  it should "successfully parse a main.wdl that imports the protected doc when the token is correct" in {
    val (doc, _, _) = parseMain(resolverWith(Map(serverHost -> ExpectedToken)))
    doc.workflow.map(_.name) shouldBe Some("main")
  }
}
