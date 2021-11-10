package dx.core

import dx.api.ExecutionEnvironment
import dx.core.ir.{DxName, SimpleDxName}

object Constants {
  // Key used in object values to signify it is a non-native dx
  // value that needs to be unwrapped before being deserialized.
  // This sequence is also used as a suffix to denote a reserved
  // parameter name.
  val ComplexValueKey = "___"
  // Suffix used for field names that are secondary to
  // complex-valued field names. A parameter with this suffix is
  // always of type array:file and is used to enumerate the files
  // that are nested within the complex value hash.
  val FlatFilesSuffix = "___dxfiles"

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
  val ScatterChunkIndex = "scatterChunkIndex"
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
  val OutputShape = "output_shape___"
  val SkippedIndices = "skipped_indices___"

  // stages that the compiler uses in generated DNAx workflows
  val CommonStage = "common"
  val EvalStage = "eval"
  val ReorgStage = "reorg"
  val OutputStage = "outputs"

  private def parameterName(prefix: String): SimpleDxName = {
    SimpleDxName.fromSourceName(prefix, Some(ComplexValueKey))
  }

  // reserved parameter names
  val InputManifest: DxName = parameterName("input_manifest")
  val InputManifestFiles: DxName = parameterName("input_manifest_files")
  val InputLinks: DxName = parameterName("input_links")
  val WorkflowInputManifest: DxName = parameterName("workflow_input_manifest")
  val WorkflowInputManifestFiles: DxName = parameterName("workflow_input_manifest_files")
  val WorkflowInputLinks: DxName = parameterName("workflow_input_links")
  val OutputId: DxName = parameterName("output_id")
  val CallName: DxName = parameterName("call_name")
  val OutputManifest: DxName = parameterName("output_manifest")
  val ExtraOutputs: DxName = parameterName("extra_outputs")
  val Overrides: DxName = parameterName("overrides")
  val ValueKey: DxName = parameterName("value")
  val WorkflowKey: DxName = parameterName("workflow")
  val ReorgConfig: DxName = parameterName("reorg_conf")
  val ReorgStatus: DxName = parameterName("reorg_status")
  // value of reorg_status output when job is completed
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
  val JobsPerScatterLimit = 2000
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
