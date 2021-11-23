#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "add VEP annotation to report"
requirements:
- class: DockerRequirement
  dockerPull: "griffithlab/vatools:4.1.0"
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
      position: 1
  vep_fields:
    type:
      type: array
      items: string
    default: ['Consequence', 'SYMBOL', 'Feature', 'HGVSc', 'HGVSp']
    inputBinding:
      position: 2
  tsv:
    type: File?
    inputBinding:
      prefix: "-t"
  prefix:
    type: string?
    default: 'variants'
baseCommand: ["vep-annotation-reporter"]
arguments: ["-o", valueFrom: $(runtime.outdir)/$(inputs.prefix).annotated.tsv]
outputs:
  annotated_variants_tsv:
    type: File
    outputBinding:
      glob: $(inputs.prefix).annotated.tsv
