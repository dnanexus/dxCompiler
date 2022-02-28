package dx.compiler

import dx.api.DxApi
import dx.core.ir.{InstanceTypeSelection, ParameterLinkSerializer}
import dx.util.{FileAccessProtocol, FileSourceResolver, Logger}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.mockito.MockitoSugar.mock

class WorkflowCompilerTest extends AnyFlatSpec with Matchers {
  "WorkflowCompiler" should "generate a JSON output spec" in {

    val mockUserProtocols = mock[FileAccessProtocol]

    val fileResolver = FileSourceResolver.create(userProtocols = Vector(mockUserProtocols))

    val parameterLinkSerializer = ParameterLinkSerializer(
        fileResolver = fileResolver,
        dxApi = DxApi.get
    )

    val wfCompiler = WorkflowCompiler(
        separateOutputs = true,
        extras = None,
        parameterLinkSerializer = parameterLinkSerializer,
        useManifests = false,
        complexPathValues = false,
        fileResolver = fileResolver,
        instanceTypeSelection = InstanceTypeSelection.Static,
        dxApi = DxApi.get,
        logger = Logger.get
    )
    wfCompiler shouldBe (Some(WorkflowCompiler))
  }
}
