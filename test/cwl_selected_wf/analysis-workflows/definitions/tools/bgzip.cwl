#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "bgzip VCF"
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
  file:
    type: File
    inputBinding:
      position: 1
baseCommand: ["/usr/local/bin/bgzip"]
arguments: ["-c"]
stdout: $(inputs.file.basename).gz
outputs:
  bgzipped_file:
    type: stdout
