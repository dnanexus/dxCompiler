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
  outdirMin: 1030
  ramMin: 4096
  tmpdirMin: 3
- class: dx:InputResourceRequirement
  indirMin: 7
- class: SoftwareRequirement
  packages:
  - package: bwa
    specs:
    - https://anaconda.org/bioconda/bwa
  - package: bwakit
    specs:
    - https://anaconda.org/bioconda/bwakit
  - package: grabix
    specs:
    - https://anaconda.org/bioconda/grabix
  - package: minimap2
    specs:
    - https://anaconda.org/bioconda/minimap2
  - package: novoalign
    specs:
    - https://anaconda.org/bioconda/novoalign
  - package: snap-aligner
    specs:
    - https://anaconda.org/bioconda/snap-aligner
    version:
    - 1.0dev.97
  - package: sentieon
    specs:
    - https://anaconda.org/bioconda/sentieon
  - package: samtools
    specs:
    - https://anaconda.org/bioconda/samtools
  - package: sambamba
    specs:
    - https://anaconda.org/bioconda/sambamba
  - package: fgbio
    specs:
    - https://anaconda.org/bioconda/fgbio
  - package: umis
    specs:
    - https://anaconda.org/bioconda/umis
  - package: biobambam
    specs:
    - https://anaconda.org/bioconda/biobambam
  - package: seqtk
    specs:
    - https://anaconda.org/bioconda/seqtk
  - package: samblaster
    specs:
    - https://anaconda.org/bioconda/samblaster
  - package: variantbam
    specs:
    - https://anaconda.org/bioconda/variantbam
- class: arv:APIRequirement
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
- id: process_alignment_rec
  type:
    fields:
    - name: files
      type:
      - 'null'
      - items: File
        type: array
    - name: config__algorithm__quality_format
      type:
      - string
      - 'null'
    - name: align_split
      type:
      - string
      - 'null'
    name: process_alignment_rec
    type: record
baseCommand:
- bcbio_nextgen.py
- runfn
- process_alignment
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=single-parallel
- sentinel_outputs=work_bam,align_bam,hla__fastq,work_bam_plus__disc,work_bam_plus__sr
- sentinel_inputs=alignment_rec:record,process_alignment_rec:record
outputs:
- id: work_bam
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: align_bam
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
$namespaces:
  arv: http://arvados.org/cwl#
  dx: https://www.dnanexus.com/cwl#
