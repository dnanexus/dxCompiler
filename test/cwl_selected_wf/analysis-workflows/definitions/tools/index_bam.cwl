#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "samtools index"
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
- class: ResourceRequirement
  ramMin: 4000
- class: InitialWorkDirRequirement
  listing:
  - ${ var f = inputs.bam; delete f.secondaryFiles; return f }
- class: InlineJavascriptRequirement
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  bam:
    type: File
arguments: ["/usr/local/bin/samtools", "index", "$(runtime.outdir)/$(inputs.bam.basename)",
  "$(runtime.outdir)/$(inputs.bam.basename).bai", {valueFrom: " && ", shellQuote: false},
  "cp", "$(inputs.bam.basename).bai", "$(runtime.outdir)/$(inputs.bam.nameroot).bai"]
outputs:
  indexed_bam:
    type: File
    secondaryFiles: [.bai, ^.bai]
    outputBinding:
      glob: $(inputs.bam.basename)
