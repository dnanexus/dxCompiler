package dx.executor.wdl

import java.nio.file.{Files, Path, Paths}
import Assumptions.isLoggedIn
import Tags.{EdgeTest, NativeTest}
import dx.api._
import dx.core.{Constants, ir}
import dx.core.io.DxWorkerPaths
import dx.core.ir.{DxName, ParameterLinkSerializer, ParameterLinkValue, Type, TypeSerde}
import dx.core.languages.wdl.{WdlBlock, WdlBundle, WdlDxName, WdlUtils}
import dx.executor.{JobMeta, WorkflowAction, WorkflowExecutor}
import dx.util.{CodecUtils, FileSourceResolver, FileUtils, Logger, PosixPath}
import dx.util.protocols.DxFileAccessProtocol
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._
import wdlTools.eval.WdlValues.V_String
import wdlTools.eval.{Eval, WdlValueBindings, WdlValues}
import wdlTools.syntax.WdlVersion
import wdlTools.types.WdlTypes.T_String
import wdlTools.types.{WdlTypes, TypedAbstractSyntax => TAT}

import scala.collection.immutable.SeqMap

private case class WorkflowTestJobMeta(override val workerPaths: DxWorkerPaths,
                                       override val dxApi: DxApi = DxApi.get,
                                       override val logger: Logger = Logger.get,
                                       rawEnv: Map[DxName, (WdlTypes.T, WdlValues.V)],
                                       rawBlockPath: Vector[Int],
                                       rawInstanceTypeDb: InstanceTypeDB,
                                       rawSourceCode: String,
                                       pathToDxFile: Map[Path, DxFile] = Map.empty)
    extends JobMeta(workerPaths, WdlDxName, dxApi, logger) {
  override val project: DxProject = null

  override val rawJsInputs: Map[DxName, JsValue] =
    ParameterLinkSerializer().createFieldsFromMap(WdlUtils.toIR(rawEnv))

  override def uploadFiles(filesToUpload: Iterable[FileUpload]): Map[Path, DxFile] = {
    filesToUpload.map { upload =>
      pathToDxFile
        .collectFirst {
          case (path, dxFile) if upload.source.endsWith(path) => upload.source -> dxFile
        }
        .getOrElse {
          throw new Exception(
              s"${upload.source} does not match any of ${pathToDxFile.keys.mkString(",")}"
          )
        }
    }.toMap
  }

  override def writeRawJsOutputs(outputJs: Map[DxName, JsValue]): Unit = {}

  override val jobId: String = null

  override def runJobScriptFunction(name: String,
                                    successCodes: Option[Set[Int]] = Some(Set(0)),
                                    truncateLogs: Boolean = true,
                                    forwardStd: Boolean = false): Unit = ()

  override val analysis: Option[DxAnalysis] = None

  override val parentJob: Option[DxJob] = None

  override val instanceType: Option[String] = Some(WorkflowTestJobMeta.InstanceType)

  override def folder: String = "/"

  override def getJobDetail(name: String): Option[JsValue] = None

  override def getExecutableAttribute(name: String): Option[JsValue] = None

  private val executableDetails: Map[String, JsValue] = Map(
      Constants.BlockPath -> JsArray(rawBlockPath.map(JsNumber(_))),
      Constants.InstanceTypeDb -> JsString(
          CodecUtils.gzipAndBase64Encode(
              rawInstanceTypeDb.toJson.prettyPrint
          )
      ),
      Constants.SourceCode -> JsString(CodecUtils.gzipAndBase64Encode(rawSourceCode)),
      Constants.WfFragmentInputTypes -> TypeSerde.serializeSpec(
          WdlUtils
            .toIRTypeMap(rawEnv.map {
              case (dxName, (t, _)) => dxName -> t
            })
            .map {
              case (dxName, t) => dxName.encoded -> t
            }
      )
  )

  override def getExecutableDetail(name: String): Option[JsValue] = executableDetails.get(name)

  override def error(e: Throwable): Unit = {}
}

private object WorkflowTestJobMeta {
  val InstanceType: String = "mem_ssd_unicorn"
}

// This test module requires being logged in to the platform.
// It compiles WDL scripts without the runtime library.
// This tests the compiler Native mode, however, it creates
// dnanexus applets and workflows that are not runnable.
class WdlWorkflowExecutorTest extends AnyFlatSpec with Matchers {
  assume(isLoggedIn)
  private val logger = Logger.Quiet
  private val dxApi = DxApi()(logger)
  private val unicornInstance =
    DxInstanceType(
        WorkflowTestJobMeta.InstanceType,
        100,
        100,
        4,
        gpu = false,
        Vector(
            ExecutionEnvironment(Constants.OsDistribution,
                                 Constants.OsRelease,
                                 Vector(Constants.OsVersion))
        ),
        Some(DiskType.SSD),
        Some(1)
    )
  private val instanceTypeDB = InstanceTypeDB(
      Map(WorkflowTestJobMeta.InstanceType -> unicornInstance)
  )

  private def setup(): DxWorkerPaths = {
    // Create a clean temp directory for the task to use
    val jobRootDir: Path = Files.createTempDirectory("dxcompiler_applet_test")
    jobRootDir.toFile.deleteOnExit()
    val workerPaths = DxWorkerPaths(PosixPath(jobRootDir.toString))
    workerPaths.createCleanDirs()
    workerPaths
  }

  private def createWorkflowExecutor(
      workerPaths: DxWorkerPaths,
      sourcePath: Path,
      blockPath: Vector[Int] = Vector(0),
      env: Map[DxName, (WdlTypes.T, WdlValues.V)] = Map.empty
  ): WdlWorkflowExecutor = {
    val wfSourceCode = FileUtils.readFileContent(sourcePath)
    val jobMeta =
      WorkflowTestJobMeta(workerPaths, dxApi, logger, env, blockPath, instanceTypeDB, wfSourceCode)
    WdlWorkflowExecutor.create(jobMeta, separateOutputs = false)
  }

  private def createFileResolver(workerPaths: DxWorkerPaths): FileSourceResolver = {
    val dxProtocol = DxFileAccessProtocol(dxApi)
    FileSourceResolver.create(
        localDirectories = Vector(workerPaths.getWorkDir().asJavaPath),
        userProtocols = Vector(dxProtocol),
        logger = logger
    )
  }

  private def createEvaluator(workerPaths: DxWorkerPaths,
                              wdlVersion: WdlVersion,
                              fileResolver: FileSourceResolver): Eval = {
    Eval(workerPaths, Some(wdlVersion), Vector.empty, fileResolver, logger)
  }

  // Note: if the file doesn't exist, this throws a null pointer exception
  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  private def parse(path: Path): WdlBundle = {
    val (doc, _) = WdlUtils.parseAndCheckSourceFile(path)
    WdlBundle.create(doc)
  }

  it should "second block in a linear workflow" in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "wf_linear.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val block = subBlocks(1)
    val env: Map[String, WdlValues.V] = Map(
        "x" -> WdlValues.V_Int(3),
        "y" -> WdlValues.V_Int(5),
        "add" -> WdlValues.V_Call("add", SeqMap("result" -> WdlValues.V_Int(8)))
    )
    val decls: Vector[TAT.PrivateVariable] = block.elements.collect {
      case eNode: TAT.PrivateVariable => eNode
    }
    val expr: TAT.Expr = decls.head.expr
    val fileResolver = createFileResolver(workerPaths)
    val eval: Eval = createEvaluator(workerPaths, wdlBundle.version, fileResolver)
    val value: WdlValues.V = eval.applyExpr(expr, WdlValueBindings(env))
    value should be(WdlValues.V_Int(9))
  }

  it should "evaluate a scatter without a call" in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "scatter_no_call.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val block = subBlocks(0)

    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val env = Map.empty[DxName, (WdlTypes.T, WdlValues.V)]
    val results: Map[DxName, (WdlTypes.T, WdlValues.V)] = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe.evaluateWorkflowElementVariables(block.elements, env)
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    results.keys.map(_.decoded) should be(Set("names", "full_name"))
    results should be(
        Map(
            WdlDxName.fromSourceName("names") -> (WdlTypes.T_Array(WdlTypes.T_String,
                                                                   nonEmpty = false),
            WdlValues.V_Array(
                Vector(WdlValues.V_String("Michael"),
                       WdlValues.V_String("Lukas"),
                       WdlValues.V_String("Martin"),
                       WdlValues.V_String("Shelly"),
                       WdlValues.V_String("Amy"))
            )),
            WdlDxName.fromSourceName("full_name") ->
              (WdlTypes.T_Array(WdlTypes.T_String, nonEmpty = false),
              WdlValues.V_Array(
                  Vector(
                      WdlValues.V_String("Michael_Manhaim"),
                      WdlValues.V_String("Lukas_Manhaim"),
                      WdlValues.V_String("Martin_Manhaim"),
                      WdlValues.V_String("Shelly_Manhaim"),
                      WdlValues.V_String("Amy_Manhaim")
                  )
              ))
        )
    )
  }

  it should "evaluate a conditional without a call" in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "conditional_no_call.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val block = subBlocks.head
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val env = Map.empty[DxName, (WdlTypes.T, WdlValues.V)]
    val results: Map[DxName, (WdlTypes.T, WdlValues.V)] = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe.evaluateWorkflowElementVariables(block.elements, env)
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    results should be(
        Map(
            WdlDxName.fromSourceName("flag") -> (WdlTypes.T_Boolean,
            WdlValues.V_Boolean(true)),
            WdlDxName.fromSourceName("cats") -> (WdlTypes.T_Optional(WdlTypes.T_String),
            WdlValues.V_Optional(WdlValues.V_String("Mr. Baggins")))
        )
    )
  }

  it should "evaluate a nested conditional/scatter without a call" in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "nested_no_call.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val block = subBlocks.head
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val results: Map[DxName, (WdlTypes.T, WdlValues.V)] = wfExecutor match {
      case exe: WdlWorkflowExecutor =>
        exe.evaluateWorkflowElementVariables(block.elements,
                                             Map.empty[DxName, (WdlTypes.T, WdlValues.V)])
      case _ => throw new Exception("expected WdlWorkflowExecutor")
    }
    results should be(
        Map(
            WdlDxName.fromSourceName("z") -> (WdlTypes
              .T_Optional(WdlTypes.T_Array(WdlTypes.T_Int, nonEmpty = true)),
            WdlValues.V_Null)
        )
    )
  }

  it should "correctly reference hash input from previous frag" in {
    val path = pathFromBasename("frag_runner", "apps_1422_deref_out.wdl")

    val workerPaths = setup()
    val wfExecutor = createWorkflowExecutor(workerPaths, path, Vector(1))
    val wdlWorkflowSupport = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val wdlExpression = WdlValues.V_Array(
        Vector(
            WdlValues.V_Pair(WdlValues.V_String("Hello"), WdlValues.V_String("World")),
            WdlValues.V_Pair(WdlValues.V_String("Halo"), WdlValues.V_String("Welt"))
        )
    )
    val irX = WdlUtils.toIRValue(wdlExpression)
    val blockContext =
      wdlWorkflowSupport.evaluateBlockInputs(
          Map(WdlDxName.fromSourceName("pairs") -> (Type.TArray(Type.THash), irX))
      )
    blockContext.prereqEnv.get(WdlDxName.fromDecodedName("xx")) shouldBe Some(
        (T_String, V_String("Hello"))
    )
  }

  it should "create proper names for scatter results" in {
    val path = pathFromBasename("frag_runner", "strings.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val scatters = wf.body.collect {
      case x: TAT.Scatter => x
    }
    scatters.size shouldBe 1
    val scatterNode = scatters.head
    scatterNode.identifier should be("x")

    val workerPaths = setup()
    val wfExecutor = createWorkflowExecutor(workerPaths, path, Vector(2, 0))
    val wdlWorkflowSupport = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val wdlX = WdlValues.V_Array(Vector(WdlValues.V_String("x"), WdlValues.V_String("y")))
    val irX = WdlUtils.toIRValue(wdlX)
    val blockContext =
      wdlWorkflowSupport.evaluateBlockInputs(
          Map(WdlDxName.fromSourceName("x") -> (Type.TArray(Type.TString), irX))
      )
    blockContext.getScatterName(wdlX, 0) shouldBe "[x,y]"
  }

  it should "Make sure calls cannot be handled by evalExpressions" in {
    val workerPaths = setup()
    val path = pathFromBasename("draft2", "shapes.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val wdlWorkflowSupport = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    // Make sure an exception is thrown if eval-expressions is called with
    // a wdl-call.
    assertThrows[Exception] {
      wdlWorkflowSupport.evaluateWorkflowElementVariables(
          subBlocks(1).elements,
          Map.empty[DxName, (WdlTypes.T, WdlValues.V)]
      )
    }
  }

  it should "include only specified inputs" in {
    val path = pathFromBasename("frag_runner", "apps_1341.wdl")
    val workerPaths = setup()
    val wfExecutor = createWorkflowExecutor(workerPaths, path, Vector(1))
    wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val wdlExpression = WdlValues.V_String("hello")
    val irX = WdlUtils.toIRValue(wdlExpression)
    val blockContext =
      wfExecutor.evaluateBlockInputs(
          Map(WdlDxName.fromSourceName("wf_inp_00") -> (Type.TString, irX))
      )
    blockContext.inputEnv shouldBe Map(
        WdlDxName.fromSourceName("wf_inp_00") -> (T_String, V_String("hello"))
    )
  }

  it should "evaluate expressions that define variables" in {
    val workerPaths = setup()
    val path = pathFromBasename("draft2", "conditionals3.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val results: Map[DxName, (WdlTypes.T, WdlValues.V)] = wfExecutor match {
      case exe: WdlWorkflowExecutor =>
        exe.evaluateWorkflowElementVariables(subBlocks(0).elements,
                                             Map.empty[DxName, (WdlTypes.T, WdlValues.V)])
      case _ => throw new Exception("expected WdlWorkflowExecutor")
    }
    results.keys.map(_.decoded) should be(Set("powers10", "i1", "i2", "i3"))
    results(WdlDxName.fromSourceName("i1")) should be(
        (WdlTypes.T_Optional(WdlTypes.T_Int), WdlValues.V_Optional(WdlValues.V_Int(1)))
    )
    results(WdlDxName.fromSourceName("i2")) should be(
        (WdlTypes.T_Optional(WdlTypes.T_Int), WdlValues.V_Null)
    )
    results(WdlDxName.fromSourceName("i3")) should be(
        (WdlTypes.T_Optional(WdlTypes.T_Int), WdlValues.V_Optional(WdlValues.V_Int(100)))
    )
    results(WdlDxName.fromSourceName("powers10")) should be(
        (WdlTypes.T_Array(WdlTypes.T_Optional(WdlTypes.T_Int), nonEmpty = false),
         WdlValues.V_Array(
             Vector(WdlValues.V_Optional(WdlValues.V_Int(1)),
                    WdlValues.V_Null,
                    WdlValues.V_Optional(WdlValues.V_Int(100)))
         ))
    )
  }

  // find the call by recursively searching the workflow nodes and nested blocks
  private def findCallByName(callName: String, allNodes: Vector[TAT.WorkflowElement]): TAT.Call = {
    def f(nodes: Vector[TAT.WorkflowElement]): Option[TAT.Call] = {
      nodes.foldLeft(None: Option[TAT.Call]) {
        case (Some(call), _) =>
          Some(call)
        case (None, call: TAT.Call) if call.actualName == callName =>
          Some(call)
        case (None, _: TAT.Call) =>
          None
        case (None, _: TAT.PrivateVariable) =>
          None
        case (None, cond: TAT.Conditional) =>
          f(cond.body)
        case (None, scatter: TAT.Scatter) =>
          f(scatter.body)
      }
    }
    f(allNodes) match {
      case None    => throw new Exception(s"call ${callName} not found")
      case Some(x) => x
    }
  }

  it should "evaluate call inputs properly" in {
    val workerPaths = setup()
    val path = pathFromBasename("draft2", "various_calls.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val call1 = findCallByName("MaybeInt", wf.body)
    val wfSupport = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val callInputs1: Map[DxName, (WdlTypes.T, WdlValues.V)] =
      wfSupport.WdlBlockContext.evaluateCallInputs(
          call1,
          Map(WdlDxName.fromSourceName("i") -> (WdlTypes.T_Int, WdlValues.V_Int(1)))
      )
    // We need to coerce the inputs into what the callee is expecting
    callInputs1 should be(
        Map(
            WdlDxName.fromSourceName("a") -> (WdlTypes.T_Optional(WdlTypes.T_Int),
            WdlValues.V_Optional(WdlValues.V_Int(1)))
        )
    )

    val call2 = findCallByName("ManyArgs", wf.body)
    val callInputs2: Map[DxName, (WdlTypes.T, WdlValues.V)] =
      wfSupport.WdlBlockContext.evaluateCallInputs(
          call2,
          Map(
              WdlDxName.fromSourceName("powers10") -> (WdlTypes.T_Array(WdlTypes.T_Int,
                                                                        nonEmpty = false),
              WdlValues.V_Array(Vector(WdlValues.V_Int(1), WdlValues.V_Int(10))))
          )
      )

    callInputs2 should be(
        Map(
            WdlDxName.fromSourceName("a") -> (WdlTypes.T_String,
            WdlValues.V_String("hello")),
            WdlDxName.fromSourceName("b") -> (WdlTypes.T_Array(WdlTypes.T_Int),
            WdlValues.V_Array(Vector(WdlValues.V_Int(1), WdlValues.V_Int(10))))
        )
    )
  }

  it should "evaluate call inputs properly II" in {
    val workerPaths = setup()
    val path = pathFromBasename("bugs", "opt_with_default.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val call = findCallByName("Foo", wf.body)
    val wdlWfExecutor = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val callInputs: Map[DxName, (WdlTypes.T, WdlValues.V)] =
      wdlWfExecutor.WdlBlockContext.evaluateCallInputs(
          call,
          Map(
              WdlDxName
                .fromSourceName("opt") -> (WdlTypes.T_Optional(WdlTypes.T_Int), WdlValues.V_Null),
              WdlDxName.fromSourceName("ref_size") -> (WdlTypes.T_Int, WdlValues.V_Int(1))
          )
      )
    // We need to coerce the inputs into what the callee is expecting
    callInputs should be(
        Map(
            WdlDxName
              .fromSourceName("opt") -> (WdlTypes.T_Optional(WdlTypes.T_Int), WdlValues.V_Null),
            WdlDxName.fromSourceName("non_opt") -> (WdlTypes.T_Int, WdlValues.V_Int(1))
        )
    )
  }

  it should "evaluate call constant inputs" in {
    val workerPaths = setup()
    val path = pathFromBasename("nested", "two_levels.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val zincCall = findCallByName("zincWithNoParams", wf.body)
    val wfSupport = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val args = wfSupport.WdlBlockContext.evaluateCallInputs(zincCall, Map.empty)
    args shouldBe Map.empty // ("a" -> (WdlTypes.T_Int, WdlValues.V_Int(3)))
  }

  it should "expressions with structs" in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "House.wdl")
    val wdlBundle = parse(path)
    val wf: TAT.Workflow = wdlBundle.primaryCallable match {
      case Some(wf: TAT.Workflow) => wf
      case _                      => throw new Exception("unexpected")
    }
    val subBlocks = WdlBlock.createBlocks(wf.body)
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val wfSupport = wfExecutor match {
      case exe: WdlWorkflowExecutor => exe
      case _                        => throw new Exception("expected WdlWorkflowExecutor")
    }
    val results =
      wfSupport.evaluateWorkflowElementVariables(subBlocks(0).elements,
                                                 Map.empty[DxName, (WdlTypes.T, WdlValues.V)])
    results.keys shouldBe Set(
        "a",
        "b",
        "tot_height",
        "tot_num_floors",
        "streets",
        "cities",
        "tot"
    ).map(WdlDxName.fromSourceName(_))

    results(WdlDxName.fromSourceName("tot")) shouldBe (
        WdlTypes.T_Struct("House",
                          SeqMap("height" -> WdlTypes.T_Int,
                                 "num_floors" -> WdlTypes.T_Int,
                                 "street" -> WdlTypes.T_String,
                                 "city" -> WdlTypes.T_String)),
        WdlValues.V_Struct(
            "House",
            SeqMap(
                "height" -> WdlValues.V_Int(32),
                "num_floors" -> WdlValues.V_Int(4),
                "street" -> WdlValues.V_String("Alda_Mary"),
                "city" -> WdlValues.V_String("Sunnyvale_Santa Clara")
            )
        )
    )
  }

  it should "fill in missing optionals" taggedAs NativeTest in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "missing_args.wdl")
    val env: Map[DxName, (WdlTypes.T, WdlValues.V)] = Map(
        WdlDxName.fromSourceName("x") -> (WdlTypes.T_Optional(WdlTypes.T_Int), WdlValues.V_Null),
        WdlDxName.fromSourceName("y") -> (WdlTypes.T_Int, WdlValues.V_Int(5))
    )
    val wfExecutor = createWorkflowExecutor(workerPaths, path, Vector(0), env)
    val (results, msg) = wfExecutor.apply(WorkflowAction.Run)
    msg shouldBe "success Run"
    results shouldBe Map(
        WdlDxName.fromSourceName("retval") -> ParameterLinkValue(JsNumber(5), Type.TInt)
    )
  }

  it should "evaluate expressions in correct order" taggedAs NativeTest in {
    val workerPaths = setup()
    val path = pathFromBasename("frag_runner", "scatter_variable_not_found.wdl")
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val (results, msg) = wfExecutor.apply(WorkflowAction.Run)
    msg shouldBe "success Run"
    results should contain key WdlDxName.fromSourceName("bam_lane1")
    results(WdlDxName.fromSourceName("bam_lane1")) shouldBe ParameterLinkValue(
        JsArray(JsString("1_ACGT_1.bam"), JsNull),
        Type.TArray(Type.TOptional(Type.TString))
    )
    val paramName = WdlDxName.fromSourceName("bam_lane1")
    wfExecutor.jobMeta.outputSerializer
      .createFieldsFromLink(results(WdlDxName.fromSourceName("bam_lane1")), paramName) shouldBe
      Vector(paramName -> JsObject("___" -> JsArray(JsString("1_ACGT_1.bam"), JsNull)),
             paramName.withSuffix(Constants.FlatFilesSuffix) -> JsArray())
  }

  it should "handle pair field access (left/right)" taggedAs (NativeTest, EdgeTest) in {
    val path = pathFromBasename("frag_runner", "scatter_with_eval.wdl")
    val workerPaths = setup()
    val wfExecutor = createWorkflowExecutor(workerPaths, path)
    val (results, msg) = wfExecutor.apply(WorkflowAction.Run)
    msg shouldBe "success Run"
    results should contain key WdlDxName.fromSourceName("info")
    results(WdlDxName.fromSourceName("info")) shouldBe ir.ParameterLinkValue(
        JsArray(JsString("Michael_27"),
                JsString("Lukas_9"),
                JsString("Martin_13"),
                JsString("Shelly_67"),
                JsString("Amy_2")),
        Type.TArray(Type.TString)
    )
  }

  def getComplexScatterName(items: Vector[Any],
                            maxLength: Int = WorkflowExecutor.JobNameLengthLimit): String = {
    WorkflowExecutor.getComplexScatterName(items.map(i => Some(i.toString)).iterator, maxLength)
  }

  it should "Build limited sized names" in {
    getComplexScatterName(Vector(1, 2, 3)) should be("1,2,3")
    getComplexScatterName(Vector(100, 200, 300).map(_.toString), 10) shouldBe "100,200,..."
    getComplexScatterName(Vector("A", "B", "hel", "nope"), 10) shouldBe "A,B,hel,..."
    getComplexScatterName(Vector("A", "B", "C", "D", "neverland"), 17) shouldBe "A,B,C,D,neverland"
    getComplexScatterName(Vector.empty, 4) shouldBe ""
  }
}
