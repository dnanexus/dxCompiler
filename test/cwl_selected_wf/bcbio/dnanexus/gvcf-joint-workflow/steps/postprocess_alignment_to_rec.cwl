#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc
  dockerPull: quay.io/bcbio/bcbio-vc
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 1028
  ramMin: 2048
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 0
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: align_bam
  secondaryFiles:
  - .bai
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__coverage_interval
  type:
    items:
    - 'null'
    - string
    type: array
- id: config__algorithm__variant_regions
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__variant_regions_merged
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__variant_regions_orig
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__coverage
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__coverage_merged
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__coverage_orig
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__seq2c_bed_ready
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__recalibrate
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__tools_on
  type:
    items:
      items: string
      type: array
    type: array
- id: genome_resources__rnaseq__gene_bed
  type:
    items: File
    type: array
- id: genome_resources__variation__dbsnp
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: reference__twobit
  type:
    items: File
    type: array
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type:
    items: File
    type: array
- id: description
  type:
    items: string
    type: array
- id: resources
  type:
    items: string
    type: array
baseCommand:
- bcbio_nextgen.py
- runfn
- postprocess_alignment_to_rec
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=postprocess_alignment_rec:description;resources;reference__fasta__base;config__algorithm__coverage_interval;genome_resources__rnaseq__gene_bed;reference__twobit;config__algorithm__recalibrate;config__algorithm__coverage;genome_resources__variation__dbsnp;config__algorithm__tools_on;config__algorithm__variant_regions;align_bam;config__algorithm__variant_regions_merged;config__algorithm__variant_regions_orig;config__algorithm__coverage_merged;config__algorithm__coverage_orig;config__algorithm__seq2c_bed_ready
- sentinel_inputs=align_bam:var,config__algorithm__coverage_interval:var,config__algorithm__variant_regions:var,config__algorithm__variant_regions_merged:var,config__algorithm__variant_regions_orig:var,config__algorithm__coverage:var,config__algorithm__coverage_merged:var,config__algorithm__coverage_orig:var,config__algorithm__seq2c_bed_ready:var,config__algorithm__recalibrate:var,config__algorithm__tools_on:var,genome_resources__rnaseq__gene_bed:var,genome_resources__variation__dbsnp:var,reference__twobit:var,reference__fasta__base:var,description:var,resources:var
outputs:
- id: postprocess_alignment_rec
  type:
    items:
      fields:
      - name: description
        type: string
      - name: resources
        type: string
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__coverage_interval
        type:
        - 'null'
        - string
      - name: genome_resources__rnaseq__gene_bed
        type: File
      - name: reference__twobit
        type: File
      - name: config__algorithm__recalibrate
        type:
        - 'null'
        - string
        - boolean
      - name: config__algorithm__coverage
        type:
        - File
        - 'null'
      - name: genome_resources__variation__dbsnp
        type: File
      - name: config__algorithm__tools_on
        type:
          items: string
          type: array
      - name: config__algorithm__variant_regions
        type:
        - File
        - 'null'
      - name: align_bam
        type:
        - File
        - 'null'
      - name: config__algorithm__variant_regions_merged
        type:
        - File
        - 'null'
      - name: config__algorithm__variant_regions_orig
        type:
        - File
        - 'null'
      - name: config__algorithm__coverage_merged
        type:
        - File
        - 'null'
      - name: config__algorithm__coverage_orig
        type:
        - File
        - 'null'
      - name: config__algorithm__seq2c_bed_ready
        type:
        - File
        - 'null'
      name: postprocess_alignment_rec
      type: record
    type: array
$namespaces:
  dx: https://www.dnanexus.com/cwl#
