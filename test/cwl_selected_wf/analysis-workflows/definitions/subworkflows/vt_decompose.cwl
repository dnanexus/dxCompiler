#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
label: "run vt decompose"
requirements:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/vt:0.57721--hf74b74d_1
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
    secondaryFiles: [".tbi"]
baseCommand: ["vt", "decompose"]
arguments: ["-s", "-o", valueFrom: $(runtime.outdir)/decomposed.vcf.gz]
outputs:
  decomposed_vcf:
    type: File
    outputBinding:
      glob: "decomposed.vcf.gz"
