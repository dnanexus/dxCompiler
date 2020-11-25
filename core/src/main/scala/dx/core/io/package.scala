package dx.core.io

import dx.util.Enum

object StreamFiles extends Enum {
  type StreamFiles = Value
  val All, None, PerFile = Value
}
