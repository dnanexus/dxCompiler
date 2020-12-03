package dx.core.languages.cwl

import java.math.{MathContext, RoundingMode}
import dx.api.InstanceTypeRequest
import dx.core.io.DxWorkerPaths
import dx.cwl._

case class RequirementEvaluator(requirements: Vector[Requirement],
                                env: Map[String, CwlValue],
                                workerPaths: DxWorkerPaths) {
  lazy val schemaDefs: Map[String, CwlSchema] = {
    requirements
      .collect {
        case req: SchemaDefRequirement => req.asMap
      }
      .flatten
      .toMap
  }

  lazy val (javascriptEnabled, javascriptLibrary): (Boolean, Option[String]) = {
    requirements.collect {
      case req: InlineJavascriptRequirement => req
    } match {
      case Vector()    => (false, None)
      case Vector(req) => (true, req.expressionLib)
      case _ =>
        throw new Exception("found multiple InlineJavascriptRequirements")
    }
  }

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

  lazy val evaluator: Evaluator = Evaluator(javascriptEnabled, javascriptLibrary, schemaDefs)
  lazy val runtime: Runtime = Runtime.create(
      outdir = workerPaths.getOutputFilesDir(ensureExists = true),
      tmpdir = workerPaths.getTempDir(ensureExists = true)
  )
  lazy val evaluatorContext: EvaluatorContext =
    EvaluatorContext(inputs = ObjectValue(env), runtime = runtime)

  def evaluate(value: CwlValue, cwlTypes: Vector[CwlType]): CwlValue = {
    def inner(innerValue: CwlValue, innerTypes: Vector[CwlType]): CwlValue = {
      innerValue match {
        case StringValue(s) =>
          evaluator.apply(s, innerTypes, evaluatorContext)
        case ArrayValue(items) =>
          innerTypes.iterator
            .map {
              case arrayType: CwlArray =>
                try {
                  Some(ArrayValue(items.map(inner(_, arrayType.itemTypes))))
                } catch {
                  case _: Throwable => None
                }
              case _ => None
            }
            .collectFirst {
              case Some(value) => value
            }
            .getOrElse(
                throw new Exception(
                    s"array ${items} does not evaluate to any of ${innerTypes}"
                )
            )
        case ObjectValue(members) =>
          innerTypes.iterator
            .map {
              case record: CwlRecord =>
                try {
                  Some(ObjectValue(members.map {
                    case (key, value) =>
                      key -> inner(value, record.fields(key).types)
                  }))
                } catch {
                  case _: Throwable => None
                }
              case _ => None
            }
            .collectFirst {
              case Some(value) => value
            }
            .getOrElse(
                throw new Exception(
                    s"object ${members} does not evaluate to any of ${innerTypes}"
                )
            )
        case _ => value
      }
    }
    inner(value, cwlTypes)
  }

  private val MinMathContext = new MathContext(0, RoundingMode.FLOOR)
  private val MaxMathContext = new MathContext(0, RoundingMode.CEILING)
  private val MiBtoGiB = Some(1024L)
  private val CwlNumericTypes: Vector[CwlType] = Vector(CwlInt, CwlLong, CwlFloat, CwlDouble)

  private def evaluateNumeric(value: CwlValue,
                              mc: MathContext,
                              scale: Option[Long] = None): Long = {
    (evaluate(value, CwlNumericTypes), scale) match {
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
        minDiskGB = minTmpGB.getOrElse(0) + minOutGB.getOrElse(0),
        maxDiskGB = maxTmpGB.getOrElse(0) + maxOutGB.getOrElse(0),
        minCpu = resources.coresMin.map(evaluateNumeric(_, MinMathContext)),
        maxCpu = resources.coresMin.map(evaluateNumeric(_, MaxMathContext))
    )
  }
}
