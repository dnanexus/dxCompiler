package dx.core.languages.wdl

import com.sun.net.httpserver.{
  HttpExchange,
  HttpHandler,
  HttpServer,
  HttpsConfigurator,
  HttpsServer
}
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

import java.io.FileInputStream
import java.net.InetSocketAddress
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path, Paths}
import java.security.KeyStore
import javax.net.ssl.{
  HostnameVerifier,
  HttpsURLConnection,
  KeyManagerFactory,
  SSLContext,
  SSLSession,
  SSLSocketFactory,
  TrustManagerFactory
}

/**
  * End-to-end test for authenticated WDL imports.
  *
  * Exercises the two ways imports are expected to work, plus the security
  * gate around Bearer tokens:
  *   - plain HTTP, no auth required          -> success (no token attached)
  *   - plain HTTP with a token configured    -> 401 (tokens are never sent
  *                                              over cleartext HTTP)
  *   - HTTPS with a trusted self-signed cert -> success when the correct token
  *                                              is configured; 401/403 otherwise
  *
  * The self-signed certificate used by the HTTPS server is generated at
  * runtime with `keytool` and trusted by this JVM via
  * `HttpsURLConnection.setDefaultSSLSocketFactory`, mirroring how a caller
  * would inject a truststore into the JVM running dxCompiler.
  */
class WdlImportHttpAuthIntegrationTest extends AnyFlatSpec with Matchers with BeforeAndAfterAll {

  private val ExpectedToken = "secret-token-abc"
  private val LocalHost = "127.0.0.1"

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

  /** Serves the WDL document only when the correct Bearer token is presented.
    * Returns 401 when no credentials are supplied and 403 when the supplied
    * credentials do not match. */
  private val protectedHandler: HttpHandler = (exchange: HttpExchange) => {
    val auth = Option(exchange.getRequestHeaders.getFirst("Authorization"))
    if (auth.contains(s"Bearer $ExpectedToken")) {
      serveBody(exchange)
    } else {
      val status = if (auth.isEmpty) 401 else 403
      exchange.sendResponseHeaders(status, -1)
      exchange.close()
    }
  }

  /** Serves the WDL document to anyone, regardless of credentials. */
  private val publicHandler: HttpHandler = (exchange: HttpExchange) => serveBody(exchange)

  private def serveBody(exchange: HttpExchange): Unit = {
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
  }

  private var httpServer: HttpServer = _
  private var httpPort: Int = _
  private var httpsServer: HttpsServer = _
  private var httpsPort: Int = _
  private var tmpDir: Path = _
  private var savedSocketFactory: SSLSocketFactory = _
  private var savedVerifier: HostnameVerifier = _

  override def beforeAll(): Unit = {
    // Plain HTTP server: a protected and a public endpoint.
    httpServer = HttpServer.create(new InetSocketAddress(LocalHost, 0), 0)
    httpServer.createContext("/imported.wdl", protectedHandler)
    httpServer.createContext("/public-imported.wdl", publicHandler)
    httpServer.setExecutor(null)
    httpServer.start()
    httpPort = httpServer.getAddress.getPort

    // HTTPS server behind a runtime-generated self-signed cert.
    tmpDir = Files.createTempDirectory("wdl-import-https-")
    val sslContext = selfSignedSslContext(tmpDir)
    httpsServer = HttpsServer.create(new InetSocketAddress(LocalHost, 0), 0)
    httpsServer.setHttpsConfigurator(new HttpsConfigurator(sslContext))
    httpsServer.createContext("/imported.wdl", protectedHandler)
    httpsServer.setExecutor(null)
    httpsServer.start()
    httpsPort = httpsServer.getAddress.getPort

    // Trust the self-signed cert for the duration of this test, the same way a
    // truststore would be injected into the JVM running dxCompiler.
    savedSocketFactory = HttpsURLConnection.getDefaultSSLSocketFactory
    savedVerifier = HttpsURLConnection.getDefaultHostnameVerifier
    HttpsURLConnection.setDefaultSSLSocketFactory(sslContext.getSocketFactory)
    HttpsURLConnection.setDefaultHostnameVerifier(new HostnameVerifier {
      override def verify(hostname: String, session: SSLSession): Boolean = true
    })
  }

  override def afterAll(): Unit = {
    if (savedSocketFactory != null) {
      HttpsURLConnection.setDefaultSSLSocketFactory(savedSocketFactory)
    }
    if (savedVerifier != null) {
      HttpsURLConnection.setDefaultHostnameVerifier(savedVerifier)
    }
    if (httpServer != null) httpServer.stop(0)
    if (httpsServer != null) httpsServer.stop(0)
    if (tmpDir != null) {
      Files.walk(tmpDir).sorted(java.util.Comparator.reverseOrder()).forEach(p => Files.delete(p))
    }
  }

  /** Generate a self-signed cert with `keytool`, then build an SSLContext that
    * both presents it (server side) and trusts it (client side). */
  private def selfSignedSslContext(dir: Path): SSLContext = {
    val password = "changeit".toCharArray
    val keystorePath = dir.resolve("server.p12")
    val keytool = Paths.get(System.getProperty("java.home"), "bin", "keytool").toString
    val cmd = Seq(
        keytool,
        "-genkeypair",
        "-alias",
        "bearer",
        "-keyalg",
        "RSA",
        "-keysize",
        "2048",
        "-dname",
        "CN=localhost",
        "-ext",
        s"san=ip:${LocalHost},dns:localhost",
        "-validity",
        "1",
        "-storetype",
        "PKCS12",
        "-keystore",
        keystorePath.toString,
        "-storepass",
        "changeit",
        "-keypass",
        "changeit"
    )
    val rc = scala.sys.process.Process(cmd).!
    if (rc != 0) {
      throw new RuntimeException(s"keytool failed to generate a self-signed cert (exit ${rc})")
    }

    val keyStore = KeyStore.getInstance("PKCS12")
    val in = new FileInputStream(keystorePath.toFile)
    try keyStore.load(in, password)
    finally in.close()

    val kmf = KeyManagerFactory.getInstance(KeyManagerFactory.getDefaultAlgorithm)
    kmf.init(keyStore, password)

    // The cert lives inside a PrivateKeyEntry, which PKIX does not treat as a
    // trust anchor; copy it into a dedicated truststore as a trustedCertEntry.
    val trustStore = KeyStore.getInstance(KeyStore.getDefaultType)
    trustStore.load(null, null)
    trustStore.setCertificateEntry("bearer", keyStore.getCertificate("bearer"))
    val tmf = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm)
    tmf.init(trustStore)

    val ctx = SSLContext.getInstance("TLS")
    ctx.init(kmf.getKeyManagers, tmf.getTrustManagers, null)
    ctx
  }

  private def httpUrl: String = s"http://${LocalHost}:${httpPort}/imported.wdl"
  private def publicHttpUrl: String = s"http://${LocalHost}:${httpPort}/public-imported.wdl"
  private def httpsUrl: String = s"https://${LocalHost}:${httpsPort}/imported.wdl"

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

  // --- Direct resolver layer ------------------------------------------------

  it should "succeed over plain HTTP when the resource requires no authentication" in {
    val resolver = resolverWith(Map.empty)
    val fs = resolver.resolve(publicHttpUrl)
    new String(fs.readBytes, StandardCharsets.UTF_8) shouldBe ImportedWdl
  }

  it should "fail with HTTP 401 and surface the env var hint when no token is configured" in {
    val resolver = resolverWith(Map.empty)
    val fs = resolver.resolve(httpUrl)
    val thrown = the[Exception] thrownBy fs.readBytes
    thrown.getMessage should include("HTTP 401 Unauthorized")
    thrown.getMessage should include(WdlImportHttpAuth.TokensEnvVar)
  }

  it should "not send a bearer token over HTTP even when one is configured" in {
    val resolver = resolverWith(Map(LocalHost -> ExpectedToken))
    val fs = resolver.resolve(httpUrl)
    val thrown = the[Exception] thrownBy fs.readBytes
    thrown.getMessage should include("HTTP 401 Unauthorized")
  }

  it should "succeed over HTTPS with a trusted self-signed cert and the correct token" in {
    val resolver = resolverWith(Map(LocalHost -> ExpectedToken))
    val fs = resolver.resolve(httpsUrl)
    new String(fs.readBytes, StandardCharsets.UTF_8) shouldBe ImportedWdl
  }

  it should "fail with HTTP 401 over HTTPS when no token is configured" in {
    val resolver = resolverWith(Map.empty)
    val thrown = the[Exception] thrownBy resolver.resolve(httpsUrl).readBytes
    thrown.getMessage should include("HTTP 401 Unauthorized")
  }

  it should "fail with HTTP 403 over HTTPS when the configured token is wrong" in {
    val resolver = resolverWith(Map(LocalHost -> "not-the-right-token"))
    val thrown = the[Exception] thrownBy resolver.resolve(httpsUrl).readBytes
    thrown.getMessage should include("HTTP 403 Forbidden")
  }

  it should "only attach the bearer token to hosts present in the token map" in {
    // Wrong-host entry must not cause the right host to be authenticated.
    val resolver = resolverWith(Map("some.other.host" -> ExpectedToken))
    val thrown = the[Exception] thrownBy resolver.resolve(httpsUrl).readBytes
    thrown.getMessage should include("HTTP 401 Unauthorized")
  }

  // --- WDL parser layer: a workflow that imports the protected document -----

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

  private def parseMain(resolver: FileSourceResolver, remoteUrl: String) = {
    VersionSupport.fromSource(
        StringFileNode(mainWdlImporting(remoteUrl)),
        WdlOptions.default,
        resolver,
        DxApi()(Logger.Quiet),
        Logger.Quiet
    )
  }

  private def causeChain(thrown: Throwable): String =
    Iterator
      .iterate(Option(thrown))(_.flatMap(t => Option(t.getCause)))
      .takeWhile(_.isDefined)
      .flatMap(_.map(_.getMessage))
      .mkString(" | ")

  it should "parse a main.wdl that imports a public doc over plain HTTP" in {
    val (doc, _, _) = parseMain(resolverWith(Map.empty), publicHttpUrl)
    doc.workflow.map(_.name) shouldBe Some("main")
  }

  it should "fail to parse a main.wdl that imports the protected doc over HTTPS when no token is set" in {
    val thrown = the[Exception] thrownBy parseMain(resolverWith(Map.empty), httpsUrl)
    causeChain(thrown) should (include("401") or include("Unauthorized"))
  }

  it should "fail to parse a main.wdl that imports the protected doc over HTTPS when the token is wrong" in {
    val thrown =
      the[Exception] thrownBy parseMain(resolverWith(Map(LocalHost -> "wrong")), httpsUrl)
    causeChain(thrown) should (include("403") or include("Forbidden"))
  }

  it should "successfully parse a main.wdl that imports the protected doc over HTTPS with the correct token" in {
    val (doc, _, _) = parseMain(resolverWith(Map(LocalHost -> ExpectedToken)), httpsUrl)
    doc.workflow.map(_.name) shouldBe Some("main")
  }
}
