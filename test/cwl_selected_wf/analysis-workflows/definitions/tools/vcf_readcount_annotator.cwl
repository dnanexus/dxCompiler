#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "add bam_readcount info to vcf"
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
  bam_readcount_tsv:
    type: File
    inputBinding:
      position: 2
  data_type:
    type:
    - type: enum
      symbols: ["DNA", "RNA"]
    inputBinding:
      position: 3
  variant_type:
    type:
    - "null"
    - type: enum
      symbols: ["snv", "indel", "all"]
    inputBinding:
      prefix: "-t"
  sample_name:
    type: string?
    inputBinding:
      prefix: "-s"
baseCommand: ["vcf-readcount-annotator"]
arguments: ["-o", valueFrom: $(runtime.outdir)/annotated.bam_readcount.vcf.gz]
outputs:
  annotated_bam_readcount_vcf:
    type: File
    outputBinding:
      glob: "annotated.bam_readcount.vcf.gz"
