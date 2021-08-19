package dx.core

import dx.api.ExecutionEnvironment

object Constants {
  // keys used in details of native applets, workflows
  val ExecLinkInfo = "execLinkInfo"
  val BlockPath = "blockPath"
  val WfFragmentInputTypes = "fqnDictTypes"
  val InstanceTypeDb = "instanceTypeDB"
  val StaticInstanceType = "staticInstanceType"
  val DelayWorkspaceDestruction = "delayWorkspaceDestruction"
  val RuntimeAttributes = "runtimeAttrs"
  val PathsAsObjects = "pathsAsObjects"
  val CompilerTag = "dxCompiler"
  val SourceCode = "sourceCode"
  val ParseOptions = "parseOptions"
  val Language = "language"
  val ScatterChunkSize = "scatterChunkSize"
  val Checksum = "checksum"
  val Version = "version"
  val DockerImage = "dockerImage"
  val NetworkDockerImage = "networkDockerImage"
  val DynamicDockerImage = "dynamicDockerImage"
  val UseManifests = "useManifests"
  val FileDependencies = "fileDependencies"
  val NativeAppDependencies = "nativeAppDependencies"
  val DockerRegistryCredentialsUri = "dockerRegistryCredentialsUri"

  // keys used in details of jobs of native applets
  val ContinueStart = "continue_start___"

  // stages that the compiler uses in generated DNAx workflows
  val CommonStage = "common"
  val EvalStage = "eval"
  val ReorgStage = "reorg"
  val OutputStage = "outputs"

  // reserved parameter names
  val InputManifest = "input_manifest___"
  val InputManifestFiles = "input_manifest_files___"
  val InputLinks = "input_links___"
  val WorkflowInputManifest = "workflow_input_manifest___"
  val WorkflowInputManifestFiles = "workflow_input_manifest_files___"
  val WorkflowInputLinks = "workflow_input_links___"
  val OutputId = "output_id___"
  val CallName = "call_name___"
  val OutputManifest = "output_manifest___"
  val ValueKey = "value___"
  val WorkflowKey = "workflow___"
  val ReorgConfig = "reorg_conf___"
  val ReorgStatus = "reorg_status___"
  val ReorgStatusCompleted = "completed"

  // keys used in manifest link hashes
  val WorkflowLinksKey = "workflow"
  val StageLinksKey = "stage"

  // keys used for applet resources
  val BundledDependsKey = "bundledDepends"
  val BundledDependsNameKey = "name"
  val BundledDependsIdKey = "id"
  val BundledDependsStagesKey = "stages"

  // deprecated properties that we still need to check for old applets
  val ChecksumPropertyDeprecated = "dxCompiler_checksum"
  val VersionPropertyDeprecated = "dxCompiler_version"

  // Limits imposed on native apps.
  val JobsPerScatterLimit = 1000
  val JobPerScatterDefault = 500

  /**
    * Very long strings cause problems with bash and the UI, so we set
    * a max limit of 32k characters
    */
  val StringLengthLimit: Int = 32 * 1024

  // other constants
  val OsDistribution = "Ubuntu"
  val OsRelease = "20.04"
  val OsVersion = "0"
  val DefaultExecutionEnvironment: ExecutionEnvironment = ExecutionEnvironment(
      OsDistribution,
      OsRelease,
      Vector(OsVersion)
  )
}
