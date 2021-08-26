package dx.executor.cwl

import dx.api.DxApi
import dx.core.io.DxWorkerPaths
import dx.core.ir.DxName
import dx.core.languages.cwl.CwlDxName
import dx.executor.WorkerJobMeta
import dx.util.Logger

case class CwlJobMeta(override val workerPaths: DxWorkerPaths = DxWorkerPaths.default,
                      override val dxApi: DxApi = DxApi.get,
                      override val logger: Logger = Logger.get)
    extends WorkerJobMeta(workerPaths, dxApi, logger) {
  override protected def encodedDxName(encodedName: String): DxName = {
    CwlDxName.fromEncodedParameterName(encodedName)
  }

  override protected def decodedDxName(decodedName: String): DxName = {
    CwlDxName.fromDecodedParameterName(decodedName)
  }
}
