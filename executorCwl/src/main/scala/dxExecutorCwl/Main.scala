package dxExecutorCwl

import dx.core.io.StreamFiles
import dx.executor.{BaseCli, JobMeta}
import dx.executor.cwl.{CwlTaskExecutor, CwlWorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorCwl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  streamFiles: StreamFiles.StreamFiles): CwlTaskExecutor = {
    CwlTaskExecutor.create(meta, streamFiles)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean): CwlWorkflowExecutor = {
    CwlWorkflowExecutor.create(meta, separateOutputs)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
