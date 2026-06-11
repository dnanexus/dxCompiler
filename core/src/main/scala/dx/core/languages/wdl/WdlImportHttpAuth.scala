package dx.core.languages.wdl

import dx.util.{
  AuthenticatedHttpFileAccessProtocol,
  FileUtils,
  Logger
}

import java.nio.charset.Charset

/**
  * dxCompiler-specific factory for building an [[AuthenticatedHttpFileAccessProtocol]]
  * configured from the environment.
  *
  * Only the WDL import path uses this; CWL compilation does not support
  * http/https imports of the source document, so attaching this protocol on
  * the CWL path would have no effect.
  */
object WdlImportHttpAuth {

  /**
    * Environment variable used to provide per-domain Bearer tokens for WDL
    * imports. Format: `domain:token[;domain:token]*`, e.g.
    * `raw.githubusercontent.com:<TOKEN>;example.com:<TOKEN>`.
    */
  val TokensEnvVar: String = "DXCOMPILER_WDL_IMPORT_BEARER_TOKENS"

  /**
    * Hint surfaced in HTTP 401 errors raised by the
    * [[AuthenticatedHttpFileAccessProtocol]], telling users how to supply
    * Bearer tokens for private hosts.
    */
  val UnauthorizedHint: String =
    s"Bearer tokens can be supplied via the ${TokensEnvVar} environment variable. " +
      "The value must be in the format domain:token[;domain:token]*, for example: " +
      "raw.githubusercontent.com:<YOUR_GITHUB_TOKEN>;example.com:<YOUR_TOKEN>."

  /**
    * Builds an [[AuthenticatedHttpFileAccessProtocol]] from the
    * [[TokensEnvVar]] environment variable. Returns a protocol with an
    * empty token map (and no auth headers) when the variable is unset.
    */
  def fromEnvironment(
      encoding: Charset = FileUtils.DefaultEncoding,
      logger: Logger = Logger.Quiet
  ): AuthenticatedHttpFileAccessProtocol = {
    val domainBearerTokens = sys.env.get(TokensEnvVar) match {
      case Some(value) =>
        val parsed = AuthenticatedHttpFileAccessProtocol.parseAuthTokens(value)
        if (parsed.nonEmpty) {
          logger.trace(
              s"${TokensEnvVar} found; authenticated HTTP imports enabled for domains: ${parsed.keys
                .mkString(", ")}"
          )
        }
        parsed
      case None =>
        Map.empty[String, String]
    }
    AuthenticatedHttpFileAccessProtocol(
        encoding = encoding,
        domainBearerTokens = domainBearerTokens,
        unauthorizedHint = Some(UnauthorizedHint),
        logger = logger
    )
  }
}
