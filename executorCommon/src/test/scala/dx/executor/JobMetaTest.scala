package dx.executor

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import Assumptions.isLoggedIn
import dx.api.DxApi
import dx.core.io.DxWorkerPaths
import dx.core.ir.DxNameFactory
import dx.util.Logger
/**
abstract class JobMetaTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)
  val workerPaths: DxWorkerPaths = DxWorkerPaths()
  val dxNameFactory: DxNameFactory = DxNameFactory()
  val dxApi: DxApi = DxApi.get
  val logger: Logger = Logger.get
  val workerJobMeta = WorkerJobMeta(waitOnUpload = true)
}
**/
