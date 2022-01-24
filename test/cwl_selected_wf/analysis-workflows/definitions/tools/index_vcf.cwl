#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "vcf index"
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
- class: ResourceRequirement
  ramMin: 4000
- class: StepInputExpressionRequirement
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  vcf:
    type: File
    inputBinding:
      valueFrom: $(self.basename)
      position: 1
arguments: ["cp", $(inputs.vcf.path), "$(runtime.outdir)/$(inputs.vcf.basename)",
  {valueFrom: " && ", shellQuote: false}, "/usr/local/bin/tabix", "-p", "vcf"]
outputs:
  indexed_vcf:
    type: File
    secondaryFiles: [.tbi]
    outputBinding:
      glob: $(inputs.vcf.basename)
