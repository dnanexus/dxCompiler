package dxExecutorWdl

import dx.core.io.StreamFiles
import dx.executor.{BaseCli, JobMeta}
import dx.executor.wdl.{WdlTaskExecutor, WdlWorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorWdl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  streamFiles: StreamFiles.StreamFiles): WdlTaskExecutor = {
    WdlTaskExecutor.create(meta, streamFiles)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean): WdlWorkflowExecutor = {
    WdlWorkflowExecutor.create(meta, separateOutputs)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
