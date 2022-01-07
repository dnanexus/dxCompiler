#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "SelectVariants (GATK 4.1.8.1)"
requirements:
- class: ResourceRequirement
  ramMin: 6000
  tmpdirMin: 25000
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
      prefix: "--variant"
      position: 2
    secondaryFiles: [.tbi]
  fields:
    type:
      type: array
      items: string
      inputBinding:
        prefix: "-F"
    default: ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'set']
    inputBinding:
      position: 3
  genotype_fields:
    type:
      type: array
      items: string
      inputBinding:
        prefix: "-GF"
    default: ['GT', 'AD', 'DP', 'AF']
    inputBinding:
      position: 4
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx4g", "VariantsToTable"]
arguments: ["-O", valueFrom: $(runtime.outdir)/variants.tsv]
outputs:
  variants_tsv:
    type: File
    outputBinding:
      glob: "variants.tsv"
