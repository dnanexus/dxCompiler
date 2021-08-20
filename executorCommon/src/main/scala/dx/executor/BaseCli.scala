package dx.executor

import java.nio.file.{InvalidPathException, Paths}

import dx.core.CliUtils._
import dx.core.io.{DxWorkerPaths, StreamFiles}
import dx.util.Enum

object BaseCli {
  val MaxConcurrentUploads: Int = 8
}

abstract class BaseCli {
  val jarName: String

  def createTaskExecutor(meta: JobMeta,
                         streamFiles: StreamFiles.StreamFiles,
                         waitOnUpload: Boolean,
                         checkInstanceType: Boolean): TaskExecutor

  def createWorkflowExecutor(meta: JobMeta,
                             separateOutputs: Boolean,
                             waitOnUpload: Boolean): WorkflowExecutor[_]

  private object StreamFilesOptionSpec
      extends SingleValueOptionSpec[StreamFiles.StreamFiles](choices = StreamFiles.values.toVector) {
    override def parseValue(value: String): StreamFiles.StreamFiles =
      StreamFiles.withNameIgnoreCase(value)
  }

  private val CommonOptions: Map[String, OptionSpec] = Map(
      "waitOnUpload" -> FlagOptionSpec.default
  )
  private val TaskOptions: Map[String, OptionSpec] = CommonOptions ++ Map(
      "streamFiles" -> StreamFilesOptionSpec,
      "streamAllFiles" -> FlagOptionSpec.default,
      "checkInstanceType" -> FlagOptionSpec.default
  )
  private val WorkflowOptions: Map[String, OptionSpec] = CommonOptions ++ Map(
      "separateOutputs" -> FlagOptionSpec.default
  )

  object ExecutorKind extends Enum {
    type ExecutorKind = Value
    val Task, Workflow = Value
  }

  private[executor] def dispatchCommand(args: Vector[String]): Termination = {
    if (args.size < 3) {
      return BadUsageTermination()
    }
    val kind = ExecutorKind.withNameIgnoreCase(args(0))
    val action = args(1).replaceAll("_", "")
    val rootDir =
      try {
        Paths.get(args(2))
      } catch {
        case _: InvalidPathException =>
          return BadUsageTermination(s"${args(1)} is not a valid root directory")
      }
    val options =
      try {
        val optionSpec = if (kind == ExecutorKind.Task) TaskOptions else WorkflowOptions
        parseCommandLine(args.drop(3), optionSpec)
      } catch {
        case e: OptionParseException =>
          return BadUsageTermination("Error parsing command line options", Some(e))
      }
    val logger = initLogger(options)
    val waitOnUpload = options.getFlag("waitOnUpload")
    logger.traceLimited(s"Creating JobMeta: rootDir ${rootDir}")
    val jobMeta = WorkerJobMeta(workerPaths = DxWorkerPaths(rootDir), logger = logger)
    try {
      kind match {
        case ExecutorKind.Task =>
          val streamFiles = options.getValue[StreamFiles.StreamFiles]("streamFiles") match {
            case Some(value)                               => value
            case None if options.getFlag("streamAllFiles") => StreamFiles.All
            case None                                      => StreamFiles.PerFile
          }
          val checkInstanceType = options.getFlag("checkInstanceType")
          logger.trace(
              s"""Creating TaskExecutor: streamFiles ${streamFiles}, waitOnUpload ${waitOnUpload}, 
                 |checkInstanceType ${checkInstanceType}""".stripMargin.replaceAll("\n", " ")
          )
          val taskExecutor =
            createTaskExecutor(jobMeta, streamFiles, waitOnUpload, checkInstanceType)
          if (taskExecutor.apply()) {
            Success("task executed successfully")
          } else {
            Success("task relaunched")
          }
        case ExecutorKind.Workflow =>
          val separateOutputs = options.getFlag("separateOutputs")
          val workflowAction =
            try {
              WorkflowAction.withNameIgnoreCase(action)
            } catch {
              case _: NoSuchElementException =>
                return BadUsageTermination(s"Unknown action ${args(0)}")
            }
          logger.trace(
              s"Creating WorkflowExecutor: separateOutputs ${separateOutputs}, waitOnUpload ${waitOnUpload}"
          )
          val executor = createWorkflowExecutor(jobMeta, separateOutputs, waitOnUpload)
          val (_, successMessage) = executor.apply(workflowAction)
          Success(successMessage)
        case _ => BadUsageTermination()
      }
    } catch {
      case e: Throwable =>
        jobMeta.error(e)
        Failure(s"failure executing ${kind} action '${action}'", Some(e))
    }
  }

  private val usageMessage =
    s"""|java -jar ${jarName}.jar <task|workflow> <action> <rootdir> [options]
        |
        |Options:
        |    -streamFiles [all,none,perfile] 
        |                           Whether to mount all files with dxfuse (do not use the 
        |                           download agent), to mount no files with dxfuse (only use 
        |                           download agent), or to allow streaming to be set on a
        |                           per-file basis (the default).
        |    -separateOutputs       Whether to put output files in a separate folder based on
        |                           the job name. If not specified, then all outputs go to the
        |                           parent job's output folder.
        |    -traceLevel [0,1,2]    How much debug information to write to the
        |                           job log at runtime. Zero means write the minimum,
        |                           one is the default, and two is for internal debugging.
        |    -quiet                 Do not print warnings or informational outputs
        |    -verbose               Print detailed progress reports
        |    -verboseKey <module>   Detailed information for a specific module
        |    -logFile <path>        File to use for logging output; defaults to stderr
        |    -waitOnUpload          Whether to wait for each file upload to complete.
        |""".stripMargin

  def main(args: Vector[String]): Unit = {
    terminate(dispatchCommand(args), usageMessage)
  }
}
