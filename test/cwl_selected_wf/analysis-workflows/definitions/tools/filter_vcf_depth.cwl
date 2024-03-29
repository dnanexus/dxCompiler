#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "filter variants at sites below a given sequence depth in each sample"
requirements:
- class: DockerRequirement
  dockerPull: mgibio/depth-filter:0.1.2
- class: ResourceRequirement
  ramMin: 4000
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  vcf:
    type: File
    inputBinding:
      position: -2
  minimum_depth:
    type: int
    inputBinding:
      prefix: "--minimum_depth"
      position: -3
  sample_names:
    type: string[]
    inputBinding:
      itemSeparator: ","
      position: -1
baseCommand: ["/opt/conda/bin/python3", "/usr/bin/depth_filter.py"]
arguments: [valueFrom: "$(runtime.outdir)/depth_filtered.vcf"]
outputs:
  depth_filtered_vcf:
    type: File
    outputBinding:
      glob: "depth_filtered.vcf"
