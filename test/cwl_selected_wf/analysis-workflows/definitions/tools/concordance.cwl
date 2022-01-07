#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Concordance checking between Tumor and Normal BAM"
requirements:
- class: ResourceRequirement
  coresMin: 1
  ramMin: 8000
  tmpdirMin: 10000
- class: DockerRequirement
  dockerPull: "brentp/somalier:v0.1.5"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  vcf:
    type: File
    inputBinding:
      prefix: "-s"
      position: 1
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "-f"
      position: 2
  bam_1:
    type: File
    inputBinding:
      position: 3
    secondaryFiles: [.bai]
  bam_2:
    type: File
    inputBinding:
      position: 4
    secondaryFiles: [.bai]
  bam_3:
    type: File?
    inputBinding:
      position: 5
    secondaryFiles: [.bai]
baseCommand: ["/usr/bin/somalier"]
arguments: ["-o", "concordance"]
outputs:
  somalier_pairs:
    type: File
    outputBinding:
      glob: "concordance.somalier.pairs.tsv"
  somalier_samples:
    type: File
    outputBinding:
      glob: "concordance.somalier.samples.tsv"
