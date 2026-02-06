package dx.core.io

import dx.util.{FileAccessProtocol, FileUtils, Logger}
import java.net.URI
import java.nio.charset.Charset

/**
  * HTTP file access protocol with Bearer token authentication support.
  * 
  * Reads authentication token from WDL_IMPORT_TOKEN environment variable.
  * Only sends tokens to allowed domains (configurable via WDL_IMPORT_TOKEN_DOMAINS)
  * to prevent credential leakage to untrusted servers.
  * 
  * @param token Optional Bearer token (defaults to WDL_IMPORT_TOKEN env var)
  * @param allowedDomains Set of domains to send auth token to
  * @param encoding Character encoding for file content
  * @param logger Logger for trace/debug output
  */
case class AuthenticatedHttpFileAccessProtocol(
    token: Option[String] = None,
    allowedDomains: Set[String] = AuthenticatedHttpFileAccessProtocol.defaultAllowedDomains,
    encoding: Charset = FileUtils.DefaultEncoding,
    logger: Logger = Logger.Quiet
) extends FileAccessProtocol {

  override val schemes: Vector[String] = Vector(FileUtils.HttpScheme, FileUtils.HttpsScheme)
  override val supportsDirectories: Boolean = true

  /**
    * Determines if authentication should be used for the given URI.
    * Only returns true if a token is configured AND the domain is in the allowed list.
    */
  private def shouldAuthenticate(uri: URI): Boolean = {
    token.isDefined && Option(uri.getHost).exists(host =>
      allowedDomains.exists(_.equalsIgnoreCase(host))
    )
  }

  override def resolve(address: String): AuthenticatedHttpFileSource = {
    val uri = URI.create(address)
    val useAuth = shouldAuthenticate(uri)
    if (useAuth) {
      logger.trace(s"Using authenticated HTTP for import from: ${uri.getHost}")
    }
    AuthenticatedHttpFileSource(uri, encoding, isDirectory = false, if (useAuth) token else None)(address)
  }

  override def resolveDirectory(address: String): AuthenticatedHttpFileSource = {
    val uri = URI.create(address)
    val useAuth = shouldAuthenticate(uri)
    if (useAuth) {
      logger.trace(s"Using authenticated HTTP for directory import from: ${uri.getHost}")
    }
    AuthenticatedHttpFileSource(uri, encoding, isDirectory = true, if (useAuth) token else None)(address)
  }
}

object AuthenticatedHttpFileAccessProtocol {

  /** Environment variable name for the Bearer token */
  val TokenEnvVar: String = "WDL_IMPORT_TOKEN"

  /** Environment variable name for custom allowed domains */
  val DomainsEnvVar: String = "WDL_IMPORT_TOKEN_DOMAINS"

  /** Default allowed domains that will receive the auth token */
  val defaultDomains: Set[String] = Set(
    "github.com",
    "raw.githubusercontent.com"
  )

  /**
    * Gets the allowed domains from environment variable or defaults.
    * WDL_IMPORT_TOKEN_DOMAINS should be a comma-separated list of domains.
    */
  lazy val defaultAllowedDomains: Set[String] = {
    sys.env.get(DomainsEnvVar) match {
      case Some(domains) =>
        domains.split(",").map(_.trim.toLowerCase).filter(_.nonEmpty).toSet
      case None =>
        defaultDomains
    }
  }

  /**
    * Creates an instance with configuration from environment variables.
    * 
    * @param logger Logger for trace output (token values are never logged)
    * @return AuthenticatedHttpFileAccessProtocol configured from environment
    */
  def fromEnvironment(logger: Logger = Logger.Quiet): AuthenticatedHttpFileAccessProtocol = {
    val tokenOpt = sys.env.get(TokenEnvVar)
    if (tokenOpt.isDefined) {
      logger.trace(
        s"${TokenEnvVar} found; authenticated HTTP imports enabled for domains: ${defaultAllowedDomains.mkString(", ")}"
      )
    }
    AuthenticatedHttpFileAccessProtocol(
      token = tokenOpt,
      allowedDomains = defaultAllowedDomains,
      logger = logger
    )
  }
}
