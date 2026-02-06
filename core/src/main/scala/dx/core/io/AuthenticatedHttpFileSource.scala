package dx.core.io

import dx.util.{AbstractAddressableFileNode, AddressableFileSource, FileUtils, PosixPath}
import java.io.{ByteArrayOutputStream, FileOutputStream}
import java.net.{HttpURLConnection, URI}
import java.nio.charset.Charset
import java.nio.file.{Files, Path}

/**
  * An HTTP file source that supports Bearer token authentication.
  * 
  * This class mirrors HttpFileSource from dxScala but adds support for
  * Authorization headers when accessing protected resources.
  * 
  * @param uri The URI to fetch
  * @param encoding Character encoding for reading content
  * @param isDirectory Whether this represents a directory (archive)
  * @param token Optional Bearer token for authentication
  * @param address The original address string
  */
case class AuthenticatedHttpFileSource(
    override val uri: URI,
    override val encoding: Charset,
    override val isDirectory: Boolean,
    token: Option[String]
)(override val address: String)
    extends AbstractAddressableFileNode(address, encoding) {

  private lazy val path = PosixPath(uri.getPath)

  override lazy val name: String =
    path.getName.getOrElse(throw new Exception(s"${path} is not a file"))

  override lazy val folder: String = path.getParent.map(_.toString).getOrElse("")

  override def container: String = s"${uri.getScheme}:${uri.getHost}:${folder}"

  private var hasBytes: Boolean = false

  /**
    * Execute a function with an HTTP connection, ensuring proper cleanup.
    * Adds Authorization header if a token is configured.
    */
  private def withConnection[T](fn: HttpURLConnection => T): T = {
    val url = uri.toURL
    var conn: HttpURLConnection = null
    try {
      conn = url.openConnection().asInstanceOf[HttpURLConnection]
      // Add authentication header if token is present
      token.foreach { t =>
        conn.setRequestProperty("Authorization", s"Bearer $t")
      }
      fn(conn)
    } finally {
      if (conn != null) {
        conn.disconnect()
      }
    }
  }

  override def exists: Boolean = {
    try {
      val rc = withConnection { conn =>
        conn.setRequestMethod("HEAD")
        conn.getResponseCode
      }
      rc match {
        case HttpURLConnection.HTTP_OK => true
        case HttpURLConnection.HTTP_UNAUTHORIZED =>
          throw new Exception(
            s"""HTTP 401 Unauthorized when accessing ${uri}.
               |If this is a private repository, ensure WDL_IMPORT_TOKEN is set with a valid access token.
               |For GitHub: generate a token at https://github.com/settings/tokens with 'repo' scope.
               |Current allowed domains can be configured via WDL_IMPORT_TOKEN_DOMAINS.""".stripMargin
          )
        case HttpURLConnection.HTTP_FORBIDDEN =>
          throw new Exception(
            s"""HTTP 403 Forbidden when accessing ${uri}.
               |The token may be invalid or lack the required permissions.
               |For GitHub: ensure the token has 'repo' scope for private repositories.""".stripMargin
          )
        case _ => false
      }
    } catch {
      case _: java.net.UnknownHostException => false
      case e: Exception                     => throw e
    }
  }

  override def getParent: Option[AuthenticatedHttpFileSource] = {
    if (path.getParent == null) {
      None
    } else {
      val newUri = if (isDirectory) uri.resolve("..") else uri.resolve(".")
      Some(AuthenticatedHttpFileSource(newUri, encoding, isDirectory = true, token)(newUri.toString))
    }
  }

  private def resolve(path: String, isDir: Boolean): AuthenticatedHttpFileSource = {
    val newUri = if (isDirectory) uri.resolve(path) else uri.resolve(".").resolve(path)
    AuthenticatedHttpFileSource(newUri, encoding, isDir, token)(newUri.toString)
  }

  override def resolve(path: String): AuthenticatedHttpFileSource = resolve(path, isDir = false)

  override def resolveDirectory(path: String): AuthenticatedHttpFileSource = resolve(path, isDir = true)

  override def relativize(fileSource: AddressableFileSource): String = {
    fileSource match {
      case fs: AuthenticatedHttpFileSource if isDirectory =>
        PosixPath(uri.getPath).relativize(PosixPath(fs.uri.getPath)).toString
      case fs: AuthenticatedHttpFileSource =>
        path.getParent.get.relativize(PosixPath(fs.uri.getPath)).toString
      case _ =>
        throw new Exception(s"not an AuthenticatedHttpFileSource: ${fileSource}")
    }
  }

  override lazy val size: Long = {
    try {
      withConnection(conn => conn.getContentLengthLong)
    } catch {
      case t: Throwable =>
        throw new Exception(s"Error getting size of URL ${uri}: ${t.getMessage}")
    }
  }

  private def fetchUri(buffer: java.io.OutputStream, chunkSize: Int = 16384): Int = {
    withConnection { conn =>
      val responseCode = conn.getResponseCode
      if (responseCode == HttpURLConnection.HTTP_UNAUTHORIZED) {
        throw new Exception(
          s"""HTTP 401 Unauthorized when fetching ${uri}.
             |If this is a private repository, ensure WDL_IMPORT_TOKEN is set with a valid access token.
             |For GitHub: generate a token at https://github.com/settings/tokens with 'repo' scope.""".stripMargin
        )
      } else if (responseCode == HttpURLConnection.HTTP_FORBIDDEN) {
        throw new Exception(
          s"""HTTP 403 Forbidden when fetching ${uri}.
             |The token may be invalid or lack the required permissions.""".stripMargin
        )
      } else if (responseCode != HttpURLConnection.HTTP_OK) {
        throw new Exception(s"HTTP ${responseCode} when fetching ${uri}")
      }
      
      val is = conn.getInputStream
      try {
        var nRead = 0
        var totalRead = 0
        val data = new Array[Byte](chunkSize)
        do {
          nRead = is.read(data, 0, chunkSize)
          if (nRead > 0) {
            buffer.write(data, 0, nRead)
            totalRead += nRead
          }
        } while (nRead > 0)
        totalRead
      } finally {
        is.close()
      }
    }
  }

  override lazy val readBytes: Array[Byte] = {
    checkFileSize()
    val buffer = new ByteArrayOutputStream()
    try {
      fetchUri(buffer)
      hasBytes = true
      buffer.toByteArray
    } finally {
      buffer.close()
    }
  }

  private def localizeToFile(path: Path): Unit = {
    if (hasBytes) {
      FileUtils.writeFileContent(path, new String(readBytes, encoding))
    } else {
      val buffer = new FileOutputStream(path.toFile)
      try {
        fetchUri(buffer)
      } finally {
        buffer.close()
      }
    }
  }

  override protected def localizeTo(file: Path): Unit = {
    if (isDirectory) {
      val dest = Files.createTempFile("temp", name)
      try {
        localizeToFile(dest)
        if (Files.exists(file)) {
          FileUtils.deleteRecursive(file)
        }
        FileUtils.unpackArchive(dest, file)
      } finally {
        FileUtils.deleteRecursive(dest)
      }
    } else {
      localizeToFile(file)
    }
  }

  override def isListable: Boolean = false
}
