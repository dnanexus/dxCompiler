cwlVersion: v1.0
id: bam_chrom_counter
class: Workflow
requirements:
- class: ScatterFeatureRequirement
inputs:
- id: bam
  type: File
outputs:
- id: bai
  type: File
  outputSource: slice_bam/bai
- id: count
  type: int[]
  outputSource: count_bam/count
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