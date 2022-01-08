#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "create filtered VCF"
requirements:
- class: ResourceRequirement
  ramMin: 6000
  tmpdirMin: 25000
- class: DockerRequirement
  dockerPull: "mgibio/gatk-cwl:3.6.0"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  vcf:
    type: File
    inputBinding:
      prefix: "--variant"
      position: 2
    secondaryFiles: [.tbi]
  filtered_vcf:
    type: File
    inputBinding:
      prefix: "--mask"
      position: 3
    secondaryFiles: [.tbi]
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "-R"
      position: 1
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T",
  "VariantFiltration"]
arguments: ["--maskName", "processSomatic", "--filterNotInMask", "-o", valueFrom: $(runtime.outdir)/output.vcf.gz]
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: "output.vcf.gz"
