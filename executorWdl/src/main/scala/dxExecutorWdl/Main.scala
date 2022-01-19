package dxExecutorWdl

import dx.core.io.StreamFiles
import dx.core.ir.DxNameFactory
import dx.core.languages.wdl.WdlDxName
import dx.executor.{BaseCli, JobMeta}
import dx.executor.wdl.{WdlTaskExecutor, WdlWorkflowExecutor}

object Main extends BaseCli {
  override val jarName = "dxExecutorWdl.jar"

  override protected val dxNameFactory: DxNameFactory = WdlDxName

  override def createTaskExecutor(meta: JobMeta,
                                  streamFiles: StreamFiles.StreamFiles,
                                  checkInstanceType: Boolean
  ): WdlTaskExecutor = {
    WdlTaskExecutor.create(meta, streamFiles, checkInstanceType)
  }

  override def createWorkflowExecutor(meta: JobMeta,
                                      separateOutputs: Boolean
  ): WdlWorkflowExecutor = {
    WdlWorkflowExecutor.create(meta, separateOutputs)
  }
}

object MainApp extends App {
  Main.main(args.toVector)
}
