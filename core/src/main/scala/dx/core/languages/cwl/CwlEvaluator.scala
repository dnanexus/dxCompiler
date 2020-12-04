package dx.core.languages.cwl

import dx.core.io.DxWorkerPaths
import dx.cwl._

case class CwlEvaluator(requirements: Vector[Requirement], workerPaths: DxWorkerPaths) {
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

  lazy val schemaDefs: Map[String, CwlSchema] = {
    requirements
      .collect {
        case req: SchemaDefRequirement => req.asMap
      }
      .flatten
      .toMap
  }

  lazy val runtime: Runtime = Runtime.create(
      outdir = workerPaths.getOutputFilesDir(ensureExists = true),
      tmpdir = workerPaths.getTempDir(ensureExists = true)
  )

  lazy val evaluator: Evaluator = Evaluator(javascriptEnabled, javascriptLibrary, schemaDefs)

  def createEvauatorContext(env: Map[String, (CwlType, CwlValue)] = Map.empty,
                            self: CwlValue = NullValue): EvaluatorContext = {
    val values = env.map {
      case (key, (_, value)) => key -> value
    }
    EvaluatorContext(self, ObjectValue(values), runtime)
  }

  private lazy val emptyEvaluatorContext = createEvauatorContext()

  def evaluate(value: CwlValue,
               cwlTypes: Vector[CwlType],
               ctx: EvaluatorContext = emptyEvaluatorContext): (CwlType, CwlValue) = {
    def inner(innerValue: CwlValue, innerTypes: Vector[CwlType]): (CwlType, CwlValue) = {
      innerValue match {
        case StringValue(s) =>
          evaluator.apply(s, innerTypes, ctx)
        case ArrayValue(items) =>
          innerTypes.iterator
            .map {
              case arrayType: CwlArray =>
                try {
                  val (types, values) = items.map(inner(_, arrayType.itemTypes)).unzip
                  Some((CwlArray(types.distinct), ArrayValue(values)))
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
                  val (types, values) = members.map {
                    case (key, value) =>
                      val (t, v) = inner(value, record.fields(key).types)
                      (key -> t, key -> v)
                  }.unzip
                  val recordType = CwlRecord(types.map {
                    case (name, t) => name -> CwlRecordField(name, Vector(t))
                  }.toMap)
                  Some((recordType, ObjectValue(values.toMap)))
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
        case _ => value.coerceTo(cwlTypes)
      }
    }
    inner(value, cwlTypes)
  }

  def evaluateMap(
      env: Map[String, (CwlType, CwlValue)],
      ctx: EvaluatorContext = emptyEvaluatorContext
  ): Map[String, (CwlType, CwlValue)] = {
    env.foldLeft(Map.empty[String, (CwlType, CwlValue)]) {
      case (accu, (key, (cwlType, cwlValue))) =>
        accu + (key -> evaluate(cwlValue, Vector(cwlType), ctx))
    }
  }
}
