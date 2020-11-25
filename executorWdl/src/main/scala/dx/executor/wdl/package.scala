package dx.executor

import com.typesafe.config.ConfigFactory

package object wdl {
  // The version lives in application.conf
  def getVersion: String = {
    val config = ConfigFactory.load("application.conf")
    config.getString("dxExecutorWdl.version")
  }
}
