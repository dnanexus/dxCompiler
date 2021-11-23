#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "vcf merge"
requirements:
- class: DockerRequirement
  dockerPull: mgibio/bcftools-cwl:1.12
- class: ResourceRequirement
  ramMin: 4000
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  vcfs:
    type: File[]
    inputBinding:
      position: 1
    secondaryFiles: [.tbi]
  merged_vcf_basename:
    type: string?
    default: 'merged'
baseCommand: ["/opt/bcftools/bin/bcftools", "concat"]
arguments:
- "--allow-overlaps"
- "--remove-duplicates"
- "--output-type"
- "z"
- "-o"
- {valueFrom: $(runtime.outdir)/$(inputs.merged_vcf_basename).vcf.gz}
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: $(inputs.merged_vcf_basename).vcf.gz
