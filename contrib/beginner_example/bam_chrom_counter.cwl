#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: bam_chrom_counter
requirements:
- class: ScatterFeatureRequirement
inputs:
- id: bam
  type: File
  default: "dx://project-BQbJpBj0bvygyQxgQ1800Jkk:file-FpQKQk00FgkGV3Vb3jJ8xqGV"
steps:
- id: slice_bam
  run: slice_bam.cwl
  in:
    bam: bam
  out: [bai, slices]
- id: count_bam
  run: count_bam.cwl
  scatter: bam
  in:
    bam: slice_bam/slices
  out: [count]
outputs:
- id: bai
  type: File
  outputSource: slice_bam/bai
- id: count
  type: int[]
  outputSource: count_bam/count
