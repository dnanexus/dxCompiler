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
  coresMin: 2
  outdirMin: 1035
  ramMin: 4096
  tmpdirMin: 6
- class: dx:InputResourceRequirement
  indirMin: 4
- class: SoftwareRequirement
  packages:
  - package: biobambam
    specs:
    - https://anaconda.org/bioconda/biobambam
  - package: samtools
    specs:
    - https://anaconda.org/bioconda/samtools
  - package: variantbam
    specs:
    - https://anaconda.org/bioconda/variantbam
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: alignment_rec
  type:
    fields:
    - name: description
      type: string
    - name: resources
      type: string
    - name: config__algorithm__align_split_size
      type:
      - 'null'
      - string
    - name: reference__fasta__base
      type: File
    - name: rgnames__lb
      type:
      - 'null'
      - string
    - name: rgnames__rg
      type: string
    - name: rgnames__lane
      type: string
    - name: reference__bwa__indexes
      type: File
    - name: config__algorithm__bam_clean
      type:
      - 'null'
      - string
      - boolean
    - name: files
      type:
        items: File
        type: array
    - name: config__algorithm__aligner
      type: string
    - name: rgnames__pl
      type: string
    - name: rgnames__pu
      type: string
    - name: config__algorithm__mark_duplicates
      type:
      - 'null'
      - string
      - boolean
    - name: analysis
      type: string
    - name: rgnames__sample
      type: string
    name: alignment_rec
    type: record
- id: work_bam
  secondaryFiles:
  - .bai
  type:
    items:
    - File
    - 'null'
    type: array
- id: align_bam_toolinput
  secondaryFiles:
  - .bai
  type:
    items:
    - File
    - 'null'
    type: array
- id: work_bam_plus__disc_toolinput
  secondaryFiles:
  - .bai
  type:
    items:
    - File
    - 'null'
    type: array
- id: work_bam_plus__sr_toolinput
  secondaryFiles:
  - .bai
  type:
    items:
    - File
    - 'null'
    type: array
- id: hla__fastq_toolinput
  type:
    items:
    - 'null'
    - items: File
      type: array
    type: array
baseCommand:
- bcbio_nextgen.py
- runfn
- merge_split_alignments
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=single-merge
- sentinel_outputs=align_bam,work_bam_plus__disc,work_bam_plus__sr,hla__fastq
- sentinel_inputs=alignment_rec:record,work_bam:var,align_bam:var,work_bam_plus__disc:var,work_bam_plus__sr:var,hla__fastq:var
outputs:
- id: align_bam
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: work_bam_plus__disc
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: work_bam_plus__sr
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: hla__fastq
  type:
  - 'null'
  - items: File
    type: array
$namespaces:
  dx: https://www.dnanexus.com/cwl#
