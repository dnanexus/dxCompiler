package dx.core.ir

import dx.api.DiskType.DiskType
import dx.api.{DxFile, ExecutionEnvironment, InstanceTypeRequest}
import dx.core.Constants

/**
  * Representation of the parts of dxapp.json `runSpec` that can be specified.
  */
object RunSpec {
  final case class RestartRequirement(max: Option[Long] = None,
                                      default: Option[Long] = None,
                                      errors: Map[String, Long] = Map.empty)
      extends RuntimeRequirement
  final case class TimeoutRequirement(days: Option[Long] = None,
                                      hours: Option[Long] = None,
                                      minutes: Option[Long] = None)
      extends RuntimeRequirement
  final case class IgnoreReuseRequirement(value: Boolean) extends RuntimeRequirement
  final case class AccessRequirement(network: Vector[String] = Vector.empty,
                                     project: Option[String] = None,
                                     allProjects: Option[String] = None,
                                     developer: Option[Boolean] = None,
                                     projectCreation: Option[Boolean] = None)
      extends RuntimeRequirement
  final case class StreamingRequirement(value: Boolean) extends RuntimeRequirement

  /**
    * Specification of instance type.
    *  An instance could be:
    *  Default: the platform default, useful for auxiliary calculations.
    *  Static:  instance type is known at compile time. We can start the
    *           job directly on the correct instance type.
    *  Dynamic: WDL specifies a calculation for the instance type, based
    *           on information known only at runtime. The generated app(let)
    *           will need to evalulate the expressions at runtime, and then
    *           start another job on the correct instance type.
    */
  sealed trait InstanceType
  case object DefaultInstanceType extends InstanceType
  case object DynamicInstanceType extends InstanceType
  case class StaticInstanceType(
      dxInstanceType: Option[String],
      minMemoryMB: Option[Long],
      maxMemoryMB: Option[Long],
      minDiskGB: Option[Long],
      maxDiskGB: Option[Long],
      diskType: Option[DiskType],
      minCpu: Option[Long],
      maxCpu: Option[Long],
      gpu: Option[Boolean],
      os: Option[ExecutionEnvironment]
  ) extends InstanceType {
    def toInstanceTypeRequest: InstanceTypeRequest = {
      InstanceTypeRequest(dxInstanceType,
                          minMemoryMB,
                          maxMemoryMB,
                          minDiskGB,
                          maxDiskGB,
                          diskType,
                          minCpu,
                          maxCpu,
                          gpu,
                          os.orElse(Some(Constants.DefaultExecutionEnvironment)))
    }
  }

  object StaticInstanceType {
    def apply(req: InstanceTypeRequest): StaticInstanceType = {
      StaticInstanceType(req.dxInstanceType,
                         req.minMemoryMB,
                         req.maxMemoryMB,
                         req.minDiskGB,
                         req.maxDiskGB,
                         req.diskType,
                         req.minCpu,
                         req.maxCpu,
                         req.gpu,
                         req.os)
    }
  }

  /**
    * A task may specify a container image to run under. Currently, DNAnexus only
    * supports Docker images. There are three supported options:
    *   None:    no image
    *   Network: the image resides on a network site and requires download
    *   DxFile: the image is a platform file
    * TODO: add a collection type that can contain multiple ContainerImage values
    */
  sealed trait ContainerImage
  case object NoImage extends ContainerImage
  case object NetworkDockerImage extends ContainerImage
  case class DxFileDockerImage(uri: String, tarball: DxFile) extends ContainerImage
}
