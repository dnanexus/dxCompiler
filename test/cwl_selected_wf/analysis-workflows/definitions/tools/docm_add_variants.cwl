#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
requirements:
- class: ResourceRequirement
  ramMin: 9000
  tmpdirMin: 25000
- class: DockerRequirement
  dockerPull: mgibio/gatk-cwl:3.6.0
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
  callers_vcf:
    type: File
    inputBinding:
      prefix: "--variant:callers"
      position: 2
    secondaryFiles: [.tbi]
  docm_vcf:
    type: File
    inputBinding:
      prefix: "--variant:docm"
      position: 3
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T",
  "CombineVariants"]
arguments: ["-genotypeMergeOptions", "PRIORITIZE", "--rod_priority_list", "callers,docm",
  "--setKey", "null", "-o", valueFrom: $(runtime.outdir)/merged.vcf.gz]
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: "merged.vcf.gz"
