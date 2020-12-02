package dxExecutorCwl

import dx.core.io.StreamFiles
import dx.executor.cwl.CwlTaskExecutor
import dx.executor.{BaseCli, FileUploader, JobMeta, WorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorWdl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  fileUploader: FileUploader,
                                  streamFiles: StreamFiles.StreamFiles): CwlTaskExecutor = {
    CwlTaskExecutor.create(meta, fileUploader, streamFiles)
  }

  override def createWorkflowExecutor(meta: JobMeta): WorkflowExecutor[_] = ???
}

object MainApp extends App {
  Main.main(args.toVector)
}
