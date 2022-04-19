cwlVersion: v1.0
id: count_bam
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
inputs:
- id: bam
  type: File
  inputBinding:
    position: 1
baseCommand: ["samtools", "view", "-c"]
outputs:
- id: count
  type: int
  outputBinding:
    glob: stdout
    loadContents: true
    outputEval: "$(parseInt(self[0].contents))"
stdout: stdout
hints:
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing