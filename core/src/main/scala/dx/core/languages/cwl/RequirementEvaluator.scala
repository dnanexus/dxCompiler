package dx.core.languages.cwl

import dx.api.{DxApi, InstanceTypeRequest}
import dx.core.Constants
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
import dx.core.ir.{InstanceTypeSelection, RuntimeRequirement}
import dx.cwl._

import scala.reflect.ClassTag

/**
  * Evaluates requirements and hints.
  * @param requirements Vector of Requirements in increasing priority order.
  *                     This means if two instances of the same requirement
  *                     are in the list, the latter takes presidence.
  * @param hints Vector of Hints in increasing priority order
  * @param env values to use for evaluation
  * @param workerPaths DxWorkerPaths
  * @param defaultRuntimeAttrs default values to use for attributes that
  *                            are not given as Requirements or Hints.
  * @param dxApi DxApi
  */
case class RequirementEvaluator(requirements: Vector[Requirement],
                                hints: Vector[Hint],
                                env: Map[String, (CwlType, CwlValue)],
                                workerPaths: DxWorkerPaths,
                                inputParameters: Map[String, InputParameter] = Map.empty,
                                defaultRuntimeAttrs: Map[String, (CwlType, CwlValue)] = Map.empty,
                                dxApi: DxApi = DxApi.get) {
  lazy val evaluator: Evaluator = Evaluator.create(requirements, hints)
  private lazy val runtime: Runtime = CwlUtils.createRuntime(workerPaths)
  private lazy val evaluatorContext: EvaluatorContext =
    CwlUtils.createEvaluatorContext(runtime,
                                    env,
                                    inputParameters = inputParameters,
                                    inputDir = workerPaths.getInputFilesDir())

  /**
    * Returns the last (i.e. highest priority) Requirement or Hint matching
    * the given filter.
    * @param filter the filter function to apply
    * @return Option[(hint, optional)]
    */
  def getHint(filter: Hint => Boolean): Option[(Hint, Boolean)] = {
    requirements.findLast(filter).map((_, false)).orElse(hints.findLast(filter).map((_, true)))
  }

  def getHintOfType[T <: Hint: ClassTag]: Option[(T, Boolean)] = {
    getHint {
      case _: T => true
      case _    => false
    } match {
      case Some((hint: T, b)) => Some((hint, b))
      case _                  => None
    }
  }

  /**
    * Returns all Requirements and Hints matching the given filter.
    * @param filter the filter function to apply
    * @return Vector of hints
    */
  def getHints(filter: Hint => Boolean): Vector[Hint] = {
    requirements.filter(filter) ++ hints.filter(filter)
  }

  def getHintsOfType[T <: Hint: ClassTag]: Vector[T] = {
    getHints({
      case _: T => true
      case _    => false
    }).map {
      case hint: T => hint
    }
  }

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

  lazy val (resources: ResourceRequirement, resourcesOptional: Boolean) = {
    getHint {
      case _: ResourceRequirement => true
      case _                      => false
    } match {
      case Some((req: ResourceRequirement, optional)) =>
        (req.merge(defaultResourceRequirement), optional)
      case _ => (defaultResourceRequirement, true)
    }
  }

  private val MinRoundingMode = BigDecimal.RoundingMode.FLOOR
  private val MaxRoundingMode = BigDecimal.RoundingMode.CEILING
  private val MiBtoGiB = Some(1024L)
  private val CwlNumericTypes = CwlMulti(Vector(CwlInt, CwlLong, CwlFloat, CwlDouble))

  private def evaluateNumeric(value: CwlValue,
                              roundingMode: BigDecimal.RoundingMode.RoundingMode,
                              scale: Option[Long] = None): Long = {
    (evaluator.evaluate(value, CwlNumericTypes, evaluatorContext)._2, scale) match {
      case (num: NumericValue, Some(d)) =>
        (num.decimalValue / d).setScale(0, roundingMode).toLongExact
      case (num: NumericValue, None) =>
        num.decimalValue.setScale(0, roundingMode).toLongExact
      case _ =>
        throw new Exception(s"${value} did not evaluate to a numeric value")
    }
  }

  def parseInstanceType: InstanceTypeRequest = {
    def getDiskGB(tmpdir: Option[CwlValue],
                  outdir: Option[CwlValue],
                  roundingMode: BigDecimal.RoundingMode.RoundingMode): Option[Long] = {
      val tmpGB = tmpdir.map(evaluateNumeric(_, roundingMode, MiBtoGiB))
      val outGB = outdir.map(evaluateNumeric(_, roundingMode, MiBtoGiB))
      if (tmpGB.isDefined || outGB.isDefined) {
        Some(tmpGB.getOrElse(0L) + outGB.getOrElse(0L))
      } else {
        None
      }
    }

    InstanceTypeRequest(
        minMemoryMB = resources.ramMin.map(evaluateNumeric(_, MinRoundingMode)),
        maxMemoryMB = resources.ramMax.map(evaluateNumeric(_, MaxRoundingMode)),
        minDiskGB =
          getDiskGB(resources.tmpdirMin, resources.outdirMin, MinRoundingMode).orElse(Some(0L)),
        maxDiskGB = getDiskGB(resources.tmpdirMax, resources.outdirMax, MaxRoundingMode),
        minCpu = resources.coresMin.map(evaluateNumeric(_, MinRoundingMode)),
        maxCpu = resources.coresMax.map(evaluateNumeric(_, MaxRoundingMode)),
        os = Some(Constants.DefaultExecutionEnvironment),
        optional = resourcesOptional
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

  def translateInstanceType(
      resolution: InstanceTypeSelection.InstanceTypeSelection
  ): InstanceType = {
    safeParseInstanceType match {
      case Some(InstanceTypeRequest.empty) => DefaultInstanceType
      case Some(req: InstanceTypeRequest)
          if resolution == InstanceTypeSelection.Dynamic && req.dxInstanceType.isDefined =>
        StaticInstanceType(InstanceTypeRequest(dxInstanceType = req.dxInstanceType))
      case Some(_) if resolution == InstanceTypeSelection.Dynamic => DynamicInstanceType
      case Some(req: InstanceTypeRequest)                         => StaticInstanceType(req)
      case None                                                   => DynamicInstanceType
    }
  }

  def translateContainer: ContainerImage = {
    getHint {
      case _: DockerRequirement => true
      case _                    => false
    } match {
      case Some((req: DockerRequirement, optional)) =>
        req.loadUri match {
          case Some(uri) if uri.startsWith("dx") =>
            val dxfile = dxApi.resolveFile(uri)
            DxFileDockerImage(uri, dxfile)
          case None if req.pullName.isDefined => NetworkDockerImage(req.pullName.get)
          case None if !optional && (req.importUri.isDefined || req.dockerfile.isDefined) =>
            throw new Exception("Docker is only supported via pull or download of a dx file")
          case _ => NoImage
        }
      case _ => NoImage
    }
  }

  def translateApplicationRequirements: Vector[RuntimeRequirement] = {
    // here we only consider the requirements that need to be applied
    // at compile time - the rest are evaluated at runtime
    Vector(
        getHint {
          case _: WorkReuseRequirement => true
          case _                       => false
        },
        getHint {
          case _: NetworkAccessRequirement => true
          case _                           => false
        },
        getHint {
          case _: ToolTimeLimitRequirement => true
          case _                           => false
        },
        getHint {
          case _: SoftwareRequirement => true
          case _                      => false
        }
    ).flatten.collect {
      case (WorkReuseRequirement(BooleanValue(allow)), _) =>
        IgnoreReuseRequirement(!allow)
      case (NetworkAccessRequirement(BooleanValue(allow)), _) =>
        // TODO: optional should have some effect here
        AccessRequirement(network = if (allow) Vector("*") else Vector.empty)
      case (ToolTimeLimitRequirement(timeLimit: NumericValue), _) =>
        TimeoutRequirement(minutes =
          Some((timeLimit.decimalValue / 60).setScale(0, MaxRoundingMode).toLongExact)
        )
      case (_: SoftwareRequirement, true) =>
        throw new Exception("SoftwareRequirement is not supported")
      case (_: InplaceUpdateRequirement, true) =>
        throw new Exception("InplaceUpdateRequirement is not supported")
    }
  }

  def translatePathDependencies: Vector[PathValue] = {
    getHints {
      case _: InitialWorkDirRequirement => true
      case _                            => false
    }.flatMap {
      case InitialWorkDirRequirement(listing) =>
        listing.collect {
          case ValueInitialWorkDirEntry(p: PathValue)     => p
          case DirInitialWorkDirEntry(p: PathValue, _, _) => p
        }
    }
  }
}
