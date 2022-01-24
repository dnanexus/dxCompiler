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
  outdirMin: 1025
  ramMin: 2048
  tmpdirMin: 1
- class: dx:InputResourceRequirement
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: htslib
    specs:
    - https://anaconda.org/bioconda/htslib
  - package: bedtools
    specs:
    - https://anaconda.org/bioconda/bedtools
  - package: pythonpy
    specs:
    - https://anaconda.org/bioconda/pythonpy
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: prep_samples_rec
  type:
    fields:
    - name: description
      type: string
    - name: resources
      type: string
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__coverage
      type:
      - 'null'
      - string
      - File
    - name: config__algorithm__variant_regions
      type: File
    name: prep_samples_rec
    type: record
baseCommand:
- bcbio_nextgen.py
- runfn
- prep_samples
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=config__algorithm__variant_regions,config__algorithm__variant_regions_merged,config__algorithm__variant_regions_orig,config__algorithm__coverage,config__algorithm__coverage_merged,config__algorithm__coverage_orig,config__algorithm__seq2c_bed_ready
- sentinel_inputs=prep_samples_rec:record
outputs:
- id: config__algorithm__variant_regions
  type:
  - File
  - 'null'
- id: config__algorithm__variant_regions_merged
  type:
  - File
  - 'null'
- id: config__algorithm__variant_regions_orig
  type:
  - File
  - 'null'
- id: config__algorithm__coverage
  type:
  - File
  - 'null'
- id: config__algorithm__coverage_merged
  type:
  - File
  - 'null'
- id: config__algorithm__coverage_orig
  type:
  - File
  - 'null'
- id: config__algorithm__seq2c_bed_ready
  type:
  - File
  - 'null'
$namespaces:
  dx: https://www.dnanexus.com/cwl#
