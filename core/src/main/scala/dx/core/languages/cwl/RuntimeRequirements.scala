package dx.core.languages.cwl

import java.math.{MathContext, RoundingMode}

import dx.api.InstanceTypeRequest
import dx.cwl._
import dx.util.Bindings

case class RuntimeRequirements(requirements: Vector[Requirement],
                               defaultAttrs: Bindings[String, CwlValue],
                               env: Map[String, CwlValue]) {

  private lazy val resourceRequirement: ResourceRequirement = {
    requirements.collect {
      case req: ResourceRequirement => req
    } match {
      case Vector()    => ResourceRequirement.default
      case Vector(req) => req.complete
      case _ =>
        throw new Exception("found multiple ResourceRequirements")
    }
  }

  private lazy val javascriptRequirement = requirements.collect {
    case req: InlineJavascriptRequirement => req
  } match {
    case Vector()    => None
    case Vector(req) => Some(req)
    case _ =>
      throw new Exception("found multiple InlineJavascriptRequirements")
  }

  private lazy val schemaDefs: Map[String, CwlSchema] = {
    requirements
      .collect {
        case req: SchemaDefRequirement => req.asMap
      }
      .flatten
      .toMap
  }

  private lazy val evaluator: Evaluator = {
    Evaluator(javascriptRequirement.isDefined,
              javascriptRequirement.flatMap(_.expressionLib),
              schemaDefs)
  }

  private val MinMathContext = new MathContext(0, RoundingMode.FLOOR)
  private val MaxMathContext = new MathContext(0, RoundingMode.CEILING)

  private def evaluateNumeric(value: CwlValue,
                              mc: MathContext,
                              scale: Option[Long] = None): Long = {
    (evaluate(value), scale) match {
      case (num: NumericValue, Some(d)) =>
        (num.decimalValue / d).round(mc).toLongExact
      case (num: NumericValue, None) =>
        num.decimalValue.round(mc).toLongExact
      case _ =>
        throw new Exception(s"${value} did not evaluate to a numeric value")
    }
  }

  def parseInstanceType: InstanceTypeRequest = {
    val minTempMb = resourceRequirement.temp
    InstanceTypeRequest(
        minMemoryMB = resourceRequirement.ramMin.map(evalaute(_, MinMathContext)),
        maxMemoryMB = resourceRequirement.ramMax.map(evalaute(_, MaxMathContext)),
      minDiskGB =
    )

  }
}
