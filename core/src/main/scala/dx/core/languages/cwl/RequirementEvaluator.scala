package dx.core.languages.cwl

import java.math.{MathContext, RoundingMode}
import dx.api.{DxApi, InstanceTypeRequest}
import dx.core.io.DxWorkerPaths
import dx.core.ir.RunSpec.{
  AccessRequirement,
  ContainerImage,
  DefaultInstanceType,
  DxFileDockerImage,
  DynamicInstanceType,
  IgnoreReuseRequirement,
  InstanceType,
  NetworkDockerImage,
  NoImage,
  StaticInstanceType,
  TimeoutRequirement
}
import dx.core.ir.RuntimeRequirement
import dx.cwl._

case class RequirementEvaluator(requirements: Vector[Requirement],
                                env: Map[String, (CwlType, CwlValue)],
                                workerPaths: DxWorkerPaths,
                                defaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = Map.empty,
                                dxApi: DxApi = DxApi.get) {
  lazy val evaluator: Evaluator = Evaluator.create(requirements)
  private lazy val runtime: Runtime = CwlUtils.createRuntime(workerPaths)
  private lazy val evaluatorContext: EvaluatorContext = CwlUtils.createEvauatorContext(runtime, env)
  private lazy val defaultResourceRequirement: ResourceRequirement = {
    val req = ResourceRequirement(
        defaultRuntimeAttrs.get("coresMin").map(_._2),
        defaultRuntimeAttrs.get("coresMax").map(_._2),
        defaultRuntimeAttrs.get("ramMin").map(_._2),
        defaultRuntimeAttrs.get("ramMax").map(_._2),
        defaultRuntimeAttrs.get("tmpdirMin").map(_._2),
        defaultRuntimeAttrs.get("tmpdirMax").map(_._2),
        defaultRuntimeAttrs.get("outdirMin").map(_._2),
        defaultRuntimeAttrs.get("outdirMax").map(_._2)
    )
    req.merge(ResourceRequirement.default)
  }

  lazy val resources: ResourceRequirement = {
    requirements.collect {
      case req: ResourceRequirement => req
    } match {
      case Vector()    => defaultResourceRequirement
      case Vector(req) => req.merge(defaultResourceRequirement)
      case _ =>
        throw new Exception("found multiple ResourceRequirements")
    }
  }

  private val MinMathContext = new MathContext(0, RoundingMode.FLOOR)
  private val MaxMathContext = new MathContext(0, RoundingMode.CEILING)
  private val MiBtoGiB = Some(1024L)
  private val CwlNumericTypes: Vector[CwlType] = Vector(CwlInt, CwlLong, CwlFloat, CwlDouble)

  private def evaluateNumeric(value: CwlValue,
                              mc: MathContext,
                              scale: Option[Long] = None): Long = {
    (evaluator.evaluate(value, CwlNumericTypes, evaluatorContext)._2, scale) match {
      case (num: NumericValue, Some(d)) =>
        (num.decimalValue / d).round(mc).toLongExact
      case (num: NumericValue, None) =>
        num.decimalValue.round(mc).toLongExact
      case _ =>
        throw new Exception(s"${value} did not evaluate to a numeric value")
    }
  }

  def parseInstanceType: InstanceTypeRequest = {
    val minTmpGB = resources.tmpdirMin.map(evaluateNumeric(_, MinMathContext, MiBtoGiB))
    val minOutGB = resources.outdirMin.map(evaluateNumeric(_, MinMathContext, MiBtoGiB))
    val maxTmpGB = resources.tmpdirMax.map(evaluateNumeric(_, MaxMathContext, MiBtoGiB))
    val maxOutGB = resources.outdirMax.map(evaluateNumeric(_, MaxMathContext, MiBtoGiB))

    InstanceTypeRequest(
        minMemoryMB = resources.ramMin.map(evaluateNumeric(_, MinMathContext)),
        maxMemoryMB = resources.ramMax.map(evaluateNumeric(_, MaxMathContext)),
        minDiskGB = Some(minTmpGB.getOrElse(0) + minOutGB.getOrElse(0)),
        maxDiskGB = Some(maxTmpGB.getOrElse(0) + maxOutGB.getOrElse(0)),
        minCpu = resources.coresMin.map(evaluateNumeric(_, MinMathContext)),
        maxCpu = resources.coresMin.map(evaluateNumeric(_, MaxMathContext))
    )
  }

  def safeParseInstanceType: Option[InstanceTypeRequest] = {
    try {
      Some(parseInstanceType)
    } catch {
      case _: Throwable =>
        // The generated code will need to calculate the instance type at runtime
        None
    }
  }

  def translateInstanceType: InstanceType = {
    safeParseInstanceType match {
      case None                            => DynamicInstanceType
      case Some(InstanceTypeRequest.empty) => DefaultInstanceType
      case Some(req: InstanceTypeRequest)  => StaticInstanceType(req)
    }
  }

  def translateContainer: ContainerImage = {
    requirements
      .collectFirst {
        case req: DockerRequirement =>
          req.loadUri match {
            case Some(uri) if uri.startsWith("dx") =>
              val dxfile = dxApi.resolveFile(uri)
              DxFileDockerImage(uri, dxfile)
            case None if req.pullName.isDefined => NetworkDockerImage
            case None if req.importUri.isDefined || req.dockerfile.isDefined =>
              throw new Exception("Docker is only supported via pull or download of a dx file")
            case _ => NoImage
          }
      }
      .getOrElse(NoImage)
  }

  def translateApplicationRequirements: Vector[RuntimeRequirement] = {
    // here we only consider the requirements that need to be applied
    // at compile time - the rest are evaluated at runtime
    requirements.collect {
      case WorkReuseRequirement(BooleanValue(allow)) =>
        IgnoreReuseRequirement(!allow)
      case NetworkAccessRequirement(BooleanValue(allow)) =>
        AccessRequirement(network = if (allow) Vector("*") else Vector.empty)
      case ToolTimeLimitRequirement(timeLimit: NumericValue) =>
        TimeoutRequirement(minutes =
          Some((timeLimit.decimalValue / 60).round(MaxMathContext).longValue)
        )
      case _: SoftwareRequirement =>
        throw new Exception("SoftwareRequirement is not supported")
    }
  }
}
