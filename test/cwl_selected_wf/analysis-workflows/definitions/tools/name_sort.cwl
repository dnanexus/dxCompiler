#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: 'sort BAM by name'
requirements:
- class: ResourceRequirement
  ramMin: 26000
  coresMin: 8
- class: DockerRequirement
  dockerPull: "mgibio/sambamba-cwl:0.6.4"
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
baseCommand: ["/usr/bin/sambamba", "sort"]
arguments: ["-t", valueFrom: $(runtime.cores), "-m", "22G", "-n", "-o", valueFrom: "$(inputs.bam.nameroot).NameSorted.bam"]
outputs:
  name_sorted_bam:
    type: File
    outputBinding:
      glob: "$(inputs.bam.nameroot).NameSorted.bam"
