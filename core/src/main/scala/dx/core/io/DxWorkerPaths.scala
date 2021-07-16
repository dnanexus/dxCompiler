package dx.core.io

import java.nio.file.{Path, Paths}
import dx.util.{BaseEvalPaths, ExecPaths, Logger}

object DxWorkerPaths {
  // The home directory on a DNAnexus worker. This directory exists only at runtime in the cloud.
  // Beware of using it in code paths that run at compile time.
  val RootDir: Path = Paths.get("/home/dnanexus")
  val TempDir = "job_scratch_space"
  val WorkDir = "work"
  val MetaDir = "meta"
  val StdoutFile = "stdout"
  val StderrFile = "stderr"
  val CommandScript = "commandScript"
  val ReturnCode = "returnCode"
  val ContainerRunScript = "containerRunScript"
  val ContainerId = "containerId"
  val InputFilesDir = "inputs"
  val OutputFilesDir = "outputs"
  val InstanceTypeDbFile = "instance_type_db.json"
  val SourceEncodedFile = "source.wdl.uu64"
  val DxdaManifestFile = "dxdaManifest.json"
  val DxfuseManifestFile = "dxfuseManifest.json"
  val DxuaManifestFile = "dxuaManifest.json"
  val DxfuseMountDir = "mnt"
  val SerializedInputsFile = "serialized_inputs.json"
  val SerializedUriToPathFile = "serialized_uri_to_path.json"
  val SerializedLocalizerFile = "serialized_localizer.json"

  lazy val default: DxWorkerPaths = DxWorkerPaths(RootDir)
}

/**
  * Important paths on a DNAnexus worker.
  * This is used in several distinct and seemingly disjoint cases:
  *  - generating applets (compile time)
  *  - localizing/delocalizing files and evaluating expressions (task/workflow runtime)
  *  - unit tests
  * Note that if a task in a container, the container will need access to temporary
  * files created with stdlib calls like "write_lines".
  * @param rootDir the root directory - typically the user's home directory
  */
case class DxWorkerPaths(rootDir: Path) extends BaseEvalPaths with ExecPaths {
  def getRootDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir("root", rootDir, ensureExists)
  }

  def getTempDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir("temp", getRootDir(ensureExists).resolve(DxWorkerPaths.TempDir), ensureExists)
  }

  /**
    * The execution directory - used as the base dir for relative paths (e.g. for glob search).
    */
  def getWorkDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir(DxWorkerPaths.WorkDir,
                   getRootDir(ensureExists).resolve(DxWorkerPaths.WorkDir),
                   ensureExists)
  }

  def getMetaDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir(DxWorkerPaths.MetaDir,
                   getRootDir(ensureExists).resolve(DxWorkerPaths.MetaDir),
                   ensureExists)
  }

  /**
    * The file that has a copy of standard output.
    */
  def getStdoutFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.StdoutFile)
  }

  /**
    * The file that has a copy of standard error.
    */
  def getStderrFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.StderrFile)
  }

  def getCommandFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.CommandScript)
  }

  def getReturnCodeFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.ReturnCode)
  }

  def getContainerCommandFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.ContainerRunScript)
  }

  def getContainerIdFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.ContainerId)
  }

  /**
    * Running applets download files from the platform to this location.
    */
  def getInputFilesDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir(DxWorkerPaths.InputFilesDir,
                   getRootDir(ensureExists).resolve(DxWorkerPaths.InputFilesDir),
                   ensureExists)
  }

  /**
    * Running applets place output files in this location.
    */
  def getOutputFilesDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir(DxWorkerPaths.OutputFilesDir,
                   getRootDir(ensureExists).resolve(DxWorkerPaths.OutputFilesDir),
                   ensureExists)
  }

  def getDxfuseMountDir(ensureExists: Boolean = false): Path = {
    getOrCreateDir(DxWorkerPaths.DxfuseMountDir,
                   getRootDir(ensureExists).resolve(DxWorkerPaths.DxfuseMountDir),
                   ensureExists)
  }

  /**
    * Where a JSON representation of the instance data base is stored.
    */
  def getInstanceTypeDbFile(ensureParentExists: Boolean = false): Path = {
    // TODO: any reason we can't put this in meta dir?
    getRootDir(ensureParentExists).resolve(DxWorkerPaths.InstanceTypeDbFile)
  }

  /**
    * Source WDL code. We could get it from the details field, but that
    * would require an additional API call. This is a private copy.
    */
  def getSourceEncodedFile(ensureParentExists: Boolean = false): Path = {
    // TODO: any reason we can't put this in meta dir?
    getRootDir(ensureParentExists).resolve(DxWorkerPaths.SourceEncodedFile)
  }

  /**
    * Location of dx download agent (dxda) manifest. It will download all these
    * files, if the file is non empty.
    */
  def getDxdaManifestFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.DxdaManifestFile)
  }

  /**
    * Location of dxfuse manifest. It will mount all these  files, if the file
    * is non empty.
    */
  def getDxfuseManifestFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.DxfuseManifestFile)
  }

  def getDxuaManifestFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.DxuaManifestFile)
  }

  /**
    * File for storing the inputs between stages of the task runner.
    */
  def getSerializedInputsFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.SerializedInputsFile)
  }

  /**
    * File for storing the dx uri to path mappings between stages of the task runner.
    */
  def getSerializedUriToPathFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.SerializedUriToPathFile)
  }

  /**
    * File for storing the LocalizationDisambiguator state between stages of the task runner.
    */
  def getSerializedLocalizerFile(ensureParentExists: Boolean = false): Path = {
    getMetaDir(ensureParentExists).resolve(DxWorkerPaths.SerializedLocalizerFile)
  }

  // create all the directory paths, so we can start using them.
  // This is used when running tasks, but NOT when compiling.
  def createCleanDirs(): Unit = {
    Logger.get.ignore(
        Vector(
            getWorkDir(ensureExists = true),
            getMetaDir(ensureExists = true),
            getTempDir(ensureExists = true),
            getInputFilesDir(ensureExists = true),
            getOutputFilesDir(ensureExists = true),
            getDxfuseMountDir(ensureExists = true)
        )
    )
  }
}
