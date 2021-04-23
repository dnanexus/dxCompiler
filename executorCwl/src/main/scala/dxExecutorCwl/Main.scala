package dxExecutorCwl

import dx.core.io.StreamFiles
import dx.executor.cwl.{CwlTaskExecutor, CwlWorkflowExecutor}
import dx.executor.{BaseCli, FileUploader, JobMeta, WorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorCwl.jar"

  override def createTaskExecutor(meta: JobMeta,
                                  fileUploader: FileUploader,
                                  streamFiles: StreamFiles.StreamFiles,
                                  waitOnUpload: Boolean): CwlTaskExecutor = {
    CwlTaskExecutor.create(meta, fileUploader, streamFiles, waitOnUpload = waitOnUpload)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean): WorkflowExecutor[_] = {
    CwlWorkflowExecutor.create(meta, separateOutputs)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
