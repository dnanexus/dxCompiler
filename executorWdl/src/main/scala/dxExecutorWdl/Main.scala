package dxExecutorWdl

import dx.core.io.StreamFiles
import dx.executor.{BaseCli, FileUploader, JobMeta}
import dx.executor.wdl.{WdlTaskExecutor, WdlWorkflowExecutor}

object Main extends BaseCli with App {
  override val jarName = "dxExecutorWdl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  fileUploader: FileUploader,
                                  streamFiles: StreamFiles.StreamFiles): WdlTaskExecutor = {
    WdlTaskExecutor.create(meta, fileUploader, streamFiles)
  }

  override def createWorkflowExecutor(meta: JobMeta): WdlWorkflowExecutor = {
    WdlWorkflowExecutor.create(meta)
  }
}
