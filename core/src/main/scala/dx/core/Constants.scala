package dx.core

object Constants {
  // keys used in details of native applets
  val ExecLinkInfo: String = "execLinkInfo"
  val BlockPath: String = "blockPath"
  val WfFragmentInputTypes: String = "fqnDictTypes"
  val InstanceTypeDb: String = "instanceTypeDB"
  val RuntimeAttributes: String = "runtimeAttrs"
  val CompilerTag = "dxCompiler"
  val SourceCode: String = "sourceCode"
  val Language: String = "language"
  val ScatterChunkSize = "scatterChunkSize"
  val Checksum = "checksum"
  val Version = "version"
  val DockerImage = "docker-image"

  // keys used in details of jobs of native applets
  val ContinueStart = "continue_start___"

  // stages that the compiler uses in generated DNAx workflows
  val CommonStage = "common"
  val EvalStage = "eval"
  val ReorgStage = "reorg"
  val OutputStage = "outputs"

  // reserved parameter names
  val InputManifests = "input_manifest___"
  val InputLinks = "input_links___"
  val OutputManifests = "output_manifests___"
  val OutputLinks = "output_links___"
  val ReorgConfig = "reorg_conf___"
  val ReorgStatus = "reorg_status___"
  val ReorgStatusCompleted = "completed"

  // keys used in manifest link hashes
  val WorkflowLinksKey = "workflow"
  val StageLinksKey = "stage"

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
}
