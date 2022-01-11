#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: 'samtools index cram'
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
- class: ResourceRequirement
  ramMin: 4000
- class: InitialWorkDirRequirement
  listing:
  - ${ var f = inputs.cram; delete f.secondaryFiles; return f }
- class: InlineJavascriptRequirement
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  cram:
    type: File
arguments: ["/usr/local/bin/samtools", "index", "$(runtime.outdir)/$(inputs.cram.basename)",
  "$(runtime.outdir)/$(inputs.cram.basename).crai", {valueFrom: " && ", shellQuote: false},
  "cp", "$(inputs.cram.basename).crai", "$(runtime.outdir)/$(inputs.cram.nameroot).crai"]
outputs:
  indexed_cram:
    type: File
    secondaryFiles: [.crai, ^.crai]
    outputBinding:
      glob: $(inputs.cram.basename)
