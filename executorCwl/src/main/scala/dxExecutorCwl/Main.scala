package dxExecutorCwl

import dx.core.io.StreamFiles
import dx.executor.{BaseCli, JobMeta}
import dx.executor.cwl.{CwlTaskExecutor, CwlWorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorCwl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  streamFiles: StreamFiles.StreamFiles,
                                  waitOnUpload: Boolean): CwlTaskExecutor = {
    CwlTaskExecutor.create(meta, streamFiles, waitOnUpload)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean,
                                      waitOnUpload: Boolean): CwlWorkflowExecutor = {
    CwlWorkflowExecutor.create(meta, separateOutputs, waitOnUpload)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
