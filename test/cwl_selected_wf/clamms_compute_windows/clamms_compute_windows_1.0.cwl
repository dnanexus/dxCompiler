baseCommand: []
class: CommandLineTool
cwlVersion: v1.0
id: clamms_compute_windows
inputs:
  genome_fasta:
    doc: FASTA file of the reference genome (optionally gzipped, default  (default
      grch37)
    inputBinding:
      position: 2
      prefix: --genome_fasta
    type: File?
  genome_fasta_index:
    doc: FASTA index file of the reference genome (.fai file, default grch37)
    inputBinding:
      position: 3
      prefix: --genome_fasta_index
    type: File?
  insert_size:
    default: 200
    doc: Sequencing insert size (in number of base pairs)
    inputBinding:
      position: 6
      prefix: --insert_size
    type: long
  mappability_bed_file:
    doc: BED file containing mappability of the reference genome (optionally gzipped)
    inputBinding:
      position: 4
      prefix: --mappability_bed_file
    type: File
  special_regions_bed_file:
    doc: BED file of known multi-copy duplication regions (expected copy number >
      3) and problematic regions to blacklist
    inputBinding:
      position: 5
      prefix: --special_regions_bed_file
    type: File
  targets_bed_file:
    doc: BED file of target regions (e.g. exons)
    inputBinding:
      position: 1
      prefix: --targets_bed_file
    type: File
label: 'CLAMMS: Compute and annotate windows from exome targets'
outputs:
  clamms_windows_bed_file:
    doc: BED file of computed and annotated CLAMMS windows
    outputBinding:
      glob: clamms_windows_bed_file/*
    type: File
requirements:
- class: DockerRequirement
  dockerOutputDirectory: /data/out
  dockerPull: pfda2dockstore/clamms_compute_windows:16
s:author:
  class: s:Person
  s:name: Evan Maxwell