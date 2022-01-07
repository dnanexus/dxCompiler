#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "samtools flagstat"
requirements:
- class: ResourceRequirement
  ramMin: 4000
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  bam:
    type: File
    inputBinding:
      position: 1
    secondaryFiles: [^.bai]
baseCommand: ["/usr/local/bin/samtools", "flagstat"]
stdout: "$(inputs.bam.basename).flagstat"
outputs:
  flagstats:
    type: stdout
