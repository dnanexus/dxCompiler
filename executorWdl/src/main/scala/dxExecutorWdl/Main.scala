package dxExecutorWdl

import dx.core.io.{DxWorkerPaths, StreamFiles}
import dx.executor.{BaseCli, JobMeta}
import dx.executor.wdl.{WdlJobMeta, WdlTaskExecutor, WdlWorkflowExecutor}
import dx.util.Logger

object Main extends BaseCli {
  override val jarName = "dxExecutorWdl.jar"

  override protected def createJobMeta(workerPaths: DxWorkerPaths, logger: Logger): WdlJobMeta = {
    WdlJobMeta(workerPaths, logger = logger)
  }

  override def createTaskExecutor(meta: JobMeta,
                                  streamFiles: StreamFiles.StreamFiles,
                                  waitOnUpload: Boolean,
                                  checkInstanceType: Boolean): WdlTaskExecutor = {
    WdlTaskExecutor.create(meta, streamFiles, waitOnUpload, checkInstanceType)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean,
                                      waitOnUpload: Boolean): WdlWorkflowExecutor = {
    WdlWorkflowExecutor.create(meta, separateOutputs, waitOnUpload)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
