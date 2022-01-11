#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "samtools sort"
requirements:
- class: ResourceRequirement
  ramMin: 4000
  coresMin: 1
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  output_filename:
    type: string
    default: sorted.bam
  input_bam:
    type: File
    inputBinding:
      position: 50
baseCommand: ["/usr/local/bin/samtools", "sort"]
arguments:
- prefix: -o
  valueFrom: $(runtime.outdir)/$(inputs.output_filename)
- prefix: -@
  valueFrom: $(runtime.cores)
outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
