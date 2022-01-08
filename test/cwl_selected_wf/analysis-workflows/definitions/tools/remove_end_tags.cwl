#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "remove END INFO tags"
requirements:
- class: ResourceRequirement
  ramMin: 4000
- class: DockerRequirement
  dockerPull: "mgibio/bcftools-cwl:1.12"
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
    secondaryFiles: [.tbi]
baseCommand: ["/opt/bcftools/bin/bcftools", "annotate"]
arguments: ["-x", "INFO/END", "-Oz", "-o", valueFrom: $(runtime.outdir)/pindel.noend.vcf.gz]
outputs:
  processed_vcf:
    type: File
    outputBinding:
      glob: "pindel.noend.vcf.gz"
