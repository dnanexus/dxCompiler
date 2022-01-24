#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "cle_annotated_vcf_filter"
requirements:
- class: DockerRequirement
  dockerPull: "mgibio/cle:v1.3.1"
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
  filter:
    type: boolean
    inputBinding:
      prefix: "filter"
      position: 1
baseCommand: ["/usr/bin/perl", "/usr/bin/docm_and_coding_indel_selection.pl"]
arguments: [$(inputs.vcf.path), $(runtime.outdir)]
outputs:
  cle_filtered_vcf:
    type: File
    outputBinding:
      glob: "annotated_filtered.vcf"
