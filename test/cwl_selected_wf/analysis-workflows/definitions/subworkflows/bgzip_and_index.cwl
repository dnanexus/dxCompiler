#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "bgzip and index VCF"
inputs:
  vcf:
    type: File
steps:
  bgzip:
    run: ../tools/bgzip.cwl
    in:
      file: vcf
    out: [bgzipped_file]
  index:
    run: ../tools/index_vcf.cwl
    in:
      vcf: bgzip/bgzipped_file
    out: [indexed_vcf]
outputs:
  indexed_vcf:
    type: File
    outputSource: index/indexed_vcf
    secondaryFiles: [.tbi]
