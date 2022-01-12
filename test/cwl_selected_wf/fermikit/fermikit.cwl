#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: fermikit
label: fermikit
requirements:
- class: DockerRequirement
  dockerOutputDirectory: /data/out
  dockerPull: pfda2dockstore/fermikit:8
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  estimated_genome_size:
    default: 3200000000
    doc: ''
    inputBinding:
      position: 5
      prefix: --estimated_genome_size
    type: long
  output_name:
    default: genome
    doc: ''
    inputBinding:
      position: 6
      prefix: --output_name
    type: string
  read_length:
    default: 150
    doc: ''
    inputBinding:
      position: 4
      prefix: --read_length
    type: long
  reads1:
    doc: gzipped FASTQ reads
    inputBinding:
      position: 1
      prefix: --reads1
    type: File
  reads2:
    doc: gzipped FASTQ reads
    inputBinding:
      position: 2
      prefix: --reads2
    type: File?
  reference_genome:
    doc: gzipped FASTA reference genome
    inputBinding:
      position: 3
      prefix: --reference_genome
    type: File
baseCommand: []
outputs:
  filtered_vcf:
    doc: ''
    outputBinding:
      glob: filtered_vcf/*
    type: File
  logs:
    doc: ''
    outputBinding:
      glob: logs/*
    type: File
  raw_vcf:
    doc: ''
    outputBinding:
      glob: raw_vcf/*
    type: File
  sv_vcf:
    doc: ''
    outputBinding:
      glob: sv_vcf/*
    type: File
  unitigs_bam:
    doc: ''
    outputBinding:
      glob: unitigs_bam/*
    type: File
s:author:
  class: s:Person
  s:name: Mike Lin
