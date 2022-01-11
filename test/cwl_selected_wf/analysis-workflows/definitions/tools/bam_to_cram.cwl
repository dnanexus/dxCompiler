#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: 'BAM to CRAM conversion'
requirements:
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
- class: ResourceRequirement
  ramMin: 4000
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "-T"
      position: 1
  bam:
    type: File
    inputBinding:
      position: 2
baseCommand: ["/usr/local/bin/samtools", "view", "-C"]
stdout: "$(inputs.bam.nameroot).cram"
outputs:
  cram:
    type: stdout
