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
  outdirMin: 1030
  ramMin: 2048
  tmpdirMin: 3
- class: dx:InputResourceRequirement
  indirMin: 0
- class: SoftwareRequirement
  packages:
  - package: multiqc
    specs:
    - https://anaconda.org/bioconda/multiqc
  - package: multiqc-bcbio
    specs:
    - https://anaconda.org/bioconda/multiqc-bcbio
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: qcout_rec
  type:
    items:
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
    type: array
baseCommand:
- bcbio_nextgen.py
- runfn
- multiqc_summary
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=summary__multiqc
- sentinel_inputs=qcout_rec:record
outputs:
- id: summary__multiqc
  type:
    items:
    - File
    - 'null'
    type: array
$namespaces:
  dx: https://www.dnanexus.com/cwl#
