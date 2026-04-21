package dx.core.io

import dx.util.{FileAccessProtocol, FileUtils, Logger}
import java.net.URI
import java.nio.charset.Charset

/**
  * HTTP file access protocol with Bearer token authentication support.
  * 
  * Reads per-domain authentication tokens from the WDL_IMPORT_TOKENS environment
  * variable. Format is semicolon-separated domain:token pairs:
  *   raw.githubusercontent.com:<TOKEN>;gitlab.com:<TOKEN>
  * 
  * Only sends tokens to domains explicitly listed in the configuration.
  * Requests to unlisted domains proceed without authentication.
  * 
  * @param domainTokens Map of domain -> Bearer token
  * @param encoding Character encoding for file content
  * @param logger Logger for trace/debug output
  */
case class AuthenticatedHttpFileAccessProtocol(
    domainTokens: Map[String, String] = Map.empty,
    encoding: Charset = FileUtils.DefaultEncoding,
    logger: Logger = Logger.Quiet
) extends FileAccessProtocol {

  override val schemes: Vector[String] = Vector(FileUtils.HttpScheme, FileUtils.HttpsScheme)
  override val supportsDirectories: Boolean = true

  /**
    * Looks up the token for the given URI's host.
    * Returns None if no token is configured for the domain.
    */
  private def tokenForUri(uri: URI): Option[String] = {
    Option(uri.getHost).flatMap { host =>
      domainTokens.collectFirst {
        case (domain, token) if domain.equalsIgnoreCase(host) => token
      }
    }
  }

  override def resolve(address: String): AuthenticatedHttpFileSource = {
    val uri = URI.create(address)
    val token = tokenForUri(uri)
    if (token.isDefined) {
      logger.trace(s"Using authenticated HTTP for import from: ${uri.getHost}")
    }
    AuthenticatedHttpFileSource(uri, encoding, isDirectory = false, token)(address)
  }

  override def resolveDirectory(address: String): AuthenticatedHttpFileSource = {
    val uri = URI.create(address)
    val token = tokenForUri(uri)
    if (token.isDefined) {
      logger.trace(s"Using authenticated HTTP for directory import from: ${uri.getHost}")
    }
    AuthenticatedHttpFileSource(uri, encoding, isDirectory = true, token)(address)
  }
}

object AuthenticatedHttpFileAccessProtocol {

  /** Environment variable name for per-domain tokens */
  val TokensEnvVar: String = "WDL_IMPORT_TOKENS"

  /**
    * Parses the WDL_IMPORT_TOKENS environment variable value.
    * Format: domain:token[;domain:token]*
    * Splits on first colon only, so tokens containing colons are supported.
    * 
    * @param value the raw env var value
    * @return Map of lowercase domain -> token
    */
  def parseTokens(value: String): Map[String, String] = {
    value
      .split(";")
      .map(_.trim)
      .filter(_.nonEmpty)
      .flatMap { entry =>
        val idx = entry.indexOf(':')
        if (idx > 0 && idx < entry.length - 1) {
          val domain = entry.substring(0, idx).trim.toLowerCase
          val token = entry.substring(idx + 1).trim
          if (domain.nonEmpty && token.nonEmpty) Some(domain -> token) else None
        } else {
          None
        }
      }
      .toMap
  }

  /**
    * Creates an instance with configuration from environment variables.
    * 
    * @param logger Logger for trace output (token values are never logged)
    * @return AuthenticatedHttpFileAccessProtocol configured from environment
    */
  def fromEnvironment(logger: Logger = Logger.Quiet): AuthenticatedHttpFileAccessProtocol = {
    val domainTokens = sys.env.get(TokensEnvVar) match {
      case Some(value) =>
        val parsed = parseTokens(value)
        if (parsed.nonEmpty) {
          logger.trace(
            s"${TokensEnvVar} found; authenticated HTTP imports enabled for domains: ${parsed.keys.mkString(", ")}"
          )
        }
        parsed
      case None =>
        Map.empty[String, String]
    }
    AuthenticatedHttpFileAccessProtocol(
      domainTokens = domainTokens,
      logger = logger
    )
  }
}
