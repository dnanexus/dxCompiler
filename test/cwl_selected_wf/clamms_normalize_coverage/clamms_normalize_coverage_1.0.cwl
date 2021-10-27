baseCommand: []
class: CommandLineTool
cwlVersion: v1.0
id: clamms_normalize_coverage
inputs:
  bam_file:
    doc: BAM file
    inputBinding:
      position: 1
      prefix: --bam_file
    type: File
  clamms_windows_bed_file:
    doc: BED file of pre-computed CLAMMS windows
    inputBinding:
      position: 2
      prefix: --clamms_windows_bed_file
    type: File
  targets_bed_file:
    doc: BED file of target regions (e.g. exons)
    inputBinding:
      position: 3
      prefix: --targets_bed_file
    type: File
label: 'CLAMMS: Calculate normalized coverage of target regions from a BAM'
outputs:
  normalized_coverage_bed_file:
    doc: Normalized read coverage summarized over CLAMMS windows
    outputBinding:
      glob: normalized_coverage_bed_file/*
    type: File
  raw_coverage_bed_file:
    doc: Raw read coverage summarized over CLAMMS windows
    outputBinding:
      glob: raw_coverage_bed_file/*
    type: File
requirements:
- class: DockerRequirement
  dockerOutputDirectory: /data/out
  dockerPull: pfda2dockstore/clamms_normalize_coverage:12
s:author:
  class: s:Person
  s:name: Evan Maxwell
