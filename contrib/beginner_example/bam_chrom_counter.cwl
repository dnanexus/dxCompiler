#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: bam_chrom_counter
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
    run: "#slice_bam"
    in:
      bam: bam
    out: [bai, slices]
  - id: count_bam
    run: "#count_bam"
    scatter: bam
    in:
      bam: slice_bam/slices
    out: [count]
- id: slice_bam
  class: CommandLineTool
  inputs:
  - id: bam
    type: File
  - id: num_chrom
    default: 22
    type: int
  outputs:
  - id: bai
    type: File
    outputBinding:
      glob: $(inputs.bam.basename).bai
  - id: slices
    type: File[]
    outputBinding:
      glob: "slices/*.bam"
  requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
  - class: InitialWorkDirRequirement
    listing:
    - entryname: slice_bam.sh
      entry: |-
        set -ex
        samtools index $1
        mkdir slices/
        for i in `seq $2`; do
            samtools view -b $1 -o slices/$i.bam chr$i
        done
    - entry: $(inputs.bam)
  baseCommand: ["sh", "slice_bam.sh"]
  arguments:
    - position: 0
      valueFrom: $(inputs.bam.basename)
    - position: 1
      valueFrom: $(inputs.num_chrom)
  hints:
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
- id: count_bam
  class: CommandLineTool
  requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
  inputs:
  - id: bam
    type: File
    inputBinding:
      position: 1
  baseCommand: ["samtools", "view", "-c"]
  outputs:
  - id: count
    type: int
    outputBinding:
      glob: stdout
      loadContents: true
      outputEval: "$(parseInt(self[0].contents))"
  stdout: stdout
  hints:
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing