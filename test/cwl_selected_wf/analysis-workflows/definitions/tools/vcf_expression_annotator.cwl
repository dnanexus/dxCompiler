#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "add expression info to vcf"
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
  expression_file:
    type: File
    inputBinding:
      position: 2
  expression_tool:
    type: string
    inputBinding:
      position: 3
  data_type:
    type: string
    inputBinding:
      position: 4
  sample_name:
    type: string
    inputBinding:
      prefix: "-s"
baseCommand: ["vcf-expression-annotator"]
arguments: ["-o", valueFrom: $(runtime.outdir)/annotated.expression.vcf.gz]
outputs:
  annotated_expression_vcf:
    type: File
    outputBinding:
      glob: "annotated.expression.vcf.gz"
