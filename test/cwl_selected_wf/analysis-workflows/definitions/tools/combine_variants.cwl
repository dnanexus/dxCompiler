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
  mutect_vcf:
    type: File
    inputBinding:
      prefix: "--variant:mutect"
      position: 2
    secondaryFiles: [.tbi]
  varscan_vcf:
    type: File
    inputBinding:
      prefix: "--variant:varscan"
      position: 3
    secondaryFiles: [.tbi]
  strelka_vcf:
    type: File
    inputBinding:
      prefix: "--variant:strelka"
      position: 4
    secondaryFiles: [.tbi]
  pindel_vcf:
    type: File
    inputBinding:
      prefix: "--variant:pindel"
      position: 5
    secondaryFiles: [.tbi]
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T",
  "CombineVariants"]
arguments: ["-genotypeMergeOptions", "PRIORITIZE", "--rod_priority_list", "mutect,varscan,strelka,pindel",
  "-o", valueFrom: $(runtime.outdir)/combined.vcf.gz]
outputs:
  combined_vcf:
    type: File
    outputBinding:
      glob: "combined.vcf.gz"
