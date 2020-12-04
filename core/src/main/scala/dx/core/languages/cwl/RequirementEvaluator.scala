package dx.core.languages.cwl

import java.math.{MathContext, RoundingMode}
import dx.api.InstanceTypeRequest
import dx.core.io.DxWorkerPaths
import dx.cwl._

case class RequirementEvaluator(requirements: Vector[Requirement],
                                env: Map[String, (CwlType, CwlValue)],
                                workerPaths: DxWorkerPaths) {
  lazy val evaluator: CwlEvaluator = CwlEvaluator(requirements, workerPaths)
  private lazy val evaluatorContext = evaluator.createEvauatorContext(env)

  lazy val resources: ResourceRequirement = {
    requirements.collect {
      case req: ResourceRequirement => req
    } match {
      case Vector()    => ResourceRequirement.default
      case Vector(req) => req.complete
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
}
