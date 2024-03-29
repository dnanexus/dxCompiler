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
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: bcftools
    specs:
    - https://anaconda.org/bioconda/bcftools
  - package: bedtools
    specs:
    - https://anaconda.org/bioconda/bedtools
  - package: fastqc
    specs:
    - https://anaconda.org/bioconda/fastqc
  - package: goleft
    specs:
    - https://anaconda.org/bioconda/goleft
  - package: hts-nim-tools
    specs:
    - https://anaconda.org/bioconda/hts-nim-tools
  - package: mosdepth
    specs:
    - https://anaconda.org/bioconda/mosdepth
  - package: picard
    specs:
    - https://anaconda.org/bioconda/picard
  - package: pythonpy
    specs:
    - https://anaconda.org/bioconda/pythonpy
  - package: qsignature
    specs:
    - https://anaconda.org/bioconda/qsignature
  - package: qualimap
    specs:
    - https://anaconda.org/bioconda/qualimap
  - package: sambamba
    specs:
    - https://anaconda.org/bioconda/sambamba
  - package: samtools
    specs:
    - https://anaconda.org/bioconda/samtools
  - package: preseq
    specs:
    - https://anaconda.org/bioconda/preseq
  - package: peddy
    specs:
    - https://anaconda.org/bioconda/peddy
- class: arv:RuntimeConstraints
  keep_cache: 4096
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: qc_rec
  type:
    fields:
    - name: description
      type: string
    - name: resources
      type: string
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__coverage_interval
      type:
      - string
      - 'null'
    - name: genome_build
      type: string
    - name: config__algorithm__coverage
      type:
      - File
      - 'null'
    - name: config__algorithm__tools_off
      type:
        items: string
        type: array
    - name: config__algorithm__qc
      type:
        items: string
        type: array
    - name: analysis
      type: string
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
    - name: config__algorithm__coverage_merged
      type:
      - File
      - 'null'
    - name: depth__samtools__stats
      type:
      - File
      - 'null'
    - name: depth__samtools__idxstats
      type:
      - File
      - 'null'
    - name: depth__variant_regions__regions
      type:
      - File
      - 'null'
    - name: depth__variant_regions__dist
      type:
      - File
      - 'null'
    - name: depth__sv_regions__regions
      type:
      - File
      - 'null'
    - name: depth__sv_regions__dist
      type:
      - File
      - 'null'
    - name: depth__coverage__regions
      type:
      - File
      - 'null'
    - name: depth__coverage__dist
      type:
      - File
      - 'null'
    - name: depth__coverage__thresholds
      type:
      - File
      - 'null'
    - name: variants__samples
      type:
        items:
          items:
          - File
          - 'null'
          type: array
        type: array
    name: qc_rec
    type: record
baseCommand:
- bcbio_nextgen.py
- runfn
- pipeline_summary
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=qcout_rec:summary__qc;summary__metrics;description;genome_build;config__algorithm__tools_off;config__algorithm__qc;config__algorithm__tools_on
- sentinel_inputs=qc_rec:record
outputs:
- id: qcout_rec
  type:
    fields:
    - name: summary__qc
      type:
      - File
      - 'null'
    - name: summary__metrics
      type:
      - string
      - 'null'
    - name: description
      type: string
    - name: genome_build
      type: string
    - name: config__algorithm__tools_off
      type:
        items: string
        type: array
    - name: config__algorithm__qc
      type:
        items: string
        type: array
    - name: config__algorithm__tools_on
      type:
        items: string
        type: array
    name: qcout_rec
    type: record
$namespaces:
  arv: http://arvados.org/cwl#
  dx: https://www.dnanexus.com/cwl#
