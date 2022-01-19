package dxExecutorCwl

import dx.core.io.StreamFiles
import dx.core.ir.DxNameFactory
import dx.core.languages.cwl.CwlDxName
import dx.executor.{BaseCli, JobMeta}
import dx.executor.cwl.{CwlTaskExecutor, CwlWorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorCwl.jar"

  override protected val dxNameFactory: DxNameFactory = CwlDxName

  override def createTaskExecutor(meta: JobMeta,
                                  streamFiles: StreamFiles.StreamFiles,
                                  checkInstanceType: Boolean
  ): CwlTaskExecutor = {
    CwlTaskExecutor.create(meta, streamFiles, checkInstanceType)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean
  ): CwlWorkflowExecutor = {
    CwlWorkflowExecutor.create(meta, separateOutputs)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
