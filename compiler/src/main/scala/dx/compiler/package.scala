package dx

import com.typesafe.config.ConfigFactory

package object compiler {
  // The version lives in application.conf
  def getVersion: String = {
    val config = ConfigFactory.load("application.conf")
    config.getString("dxCompiler.version")
  }
}
