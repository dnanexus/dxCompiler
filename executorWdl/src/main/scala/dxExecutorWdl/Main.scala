package dxExecutorWdl

import dx.core.io.StreamFiles
import dx.executor.{BaseCli, FileUploader, JobMeta}
import dx.executor.wdl.{WdlTaskExecutor, WdlWorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorWdl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  fileUploader: FileUploader,
                                  streamFiles: StreamFiles.StreamFiles,
                                  waitOnUpload: Boolean): WdlTaskExecutor = {
    WdlTaskExecutor.create(meta, fileUploader, streamFiles, waitOnUpload = waitOnUpload)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean): WdlWorkflowExecutor = {
    WdlWorkflowExecutor.create(meta, separateOutputs)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
