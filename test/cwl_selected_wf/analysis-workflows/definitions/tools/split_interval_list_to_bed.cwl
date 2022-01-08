#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 6000
- class: DockerRequirement
  dockerPull: mgibio/cle:v1.4.2
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  interval_list:
    type: File
    inputBinding:
      prefix: "INPUT="
      separate: false
      position: 1
  scatter_count:
    type: int
    inputBinding:
      prefix: "SCATTER_COUNT="
      separate: false
      position: 2
baseCommand: ['/usr/bin/perl', '/usr/bin/split_interval_list_to_bed_helper.pl']
arguments: [valueFrom: OUTPUT=$(runtime.outdir)]
outputs:
  split_beds:
    type: File[]
    outputBinding:
      glob: "*.interval.bed"
