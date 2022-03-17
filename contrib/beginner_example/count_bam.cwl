#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: count_bam
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
hints:
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: bam
  type: File
  inputBinding:
    position: 1
baseCommand: ["samtools", "view", "-c"]
stdout: stdout
outputs:
- id: count
  type: int
  outputBinding:
    glob: stdout
    loadContents: true
    outputEval: "$(parseInt(self[0].contents))"
