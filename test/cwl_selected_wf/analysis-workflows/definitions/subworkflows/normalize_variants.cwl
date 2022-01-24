#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
label: "Normalize variants"
requirements:
- class: ResourceRequirement
  ramMin: 9000
- class: DockerRequirement
  dockerPull: "broadinstitute/gatk:4.1.8.1"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "-R"
      position: 1
  vcf:
    type: File
    inputBinding:
      prefix: "-V"
      position: 2
    secondaryFiles: [".tbi"]
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx8g", "LeftAlignAndTrimVariants"]
arguments: ["-O", valueFrom: $(runtime.outdir)/normalized.vcf.gz]
outputs:
  normalized_vcf:
    type: File
    secondaryFiles: [".tbi"]
    outputBinding:
      glob: "normalized.vcf.gz"
