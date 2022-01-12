#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: 'Run optitype on dna samples'
requirements:
- class: ResourceRequirement
  ramMin: 64000
  coresMin: 4
  tmpdirMin: 20000
- class: DockerRequirement
  dockerPull: "mgibio/immuno_tools-cwl:1.0.1"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  optitype_name:
    type: string?
    default: "optitype"
    doc: "A prefix that will be used to name all files produced by the script"
    inputBinding:
      position: 1
  cram:
    type: File
    doc: "File to be HLA-typed"
    inputBinding:
      position: 2
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai]
    doc: "Reference fasta used to make the cram"
    inputBinding:
      position: 3
baseCommand: ["/bin/bash", "/usr/bin/optitype_script.sh"]
arguments: [valueFrom: $(runtime.tmpdir), valueFrom: $(runtime.outdir)]
outputs:
  optitype_tsv:
    type: File
    outputBinding:
      glob: $(inputs.optitype_name)_result.tsv
  optitype_plot:
    type: File
    outputBinding:
      glob: $(inputs.optitype_name)_coverage_plot.pdf
