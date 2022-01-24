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
  indirMin: 0
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: config__algorithm__coverage
  type:
    items:
    - File
    - 'null'
    - string
    type: array
- id: config__algorithm__variant_regions
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
- prep_samples_to_rec
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=prep_samples_rec:description;resources;reference__fasta__base;config__algorithm__coverage;config__algorithm__variant_regions
- sentinel_inputs=config__algorithm__coverage:var,config__algorithm__variant_regions:var,reference__fasta__base:var,description:var,resources:var
outputs:
- id: prep_samples_rec
  type:
    items:
      fields:
      - name: description
        type: string
      - name: resources
        type: string
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__coverage
        type:
        - File
        - 'null'
        - string
      - name: config__algorithm__variant_regions
        type: File
      name: prep_samples_rec
      type: record
    type: array
$namespaces:
  dx: https://www.dnanexus.com/cwl#
