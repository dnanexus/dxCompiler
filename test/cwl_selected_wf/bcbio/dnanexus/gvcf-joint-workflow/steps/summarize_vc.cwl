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
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: jointvc_rec
  type:
    items:
      items:
        fields:
        - name: description
          type: string
        - name: resources
          type: string
        - name: batch_samples
          type:
          - 'null'
          - items: string
            type: array
        - name: validate__summary
          type:
          - File
          - 'null'
        - name: validate__tp
          type:
          - File
          - 'null'
        - name: validate__fp
          type:
          - File
          - 'null'
        - name: validate__fn
          type:
          - File
          - 'null'
        - name: vrn_file
          type: File
        - name: config__algorithm__validate
          type:
          - File
          - 'null'
          - string
        - name: reference__fasta__base
          type: File
        - name: config__algorithm__variantcaller
          type: string
        - name: config__algorithm__coverage_interval
          type:
          - string
          - 'null'
        - name: metadata__batch
          type: string
        - name: metadata__phenotype
          type: string
        - name: reference__twobit
          type: File
        - name: reference__snpeff__hg19
          type: File
        - name: config__algorithm__validate_regions
          type:
          - File
          - 'null'
          - string
        - name: genome_build
          type: string
        - name: genome_resources__aliases__human
          type:
          - string
          - 'null'
          - boolean
        - name: config__algorithm__tools_off
          type:
            items: string
            type: array
        - name: genome_resources__variation__dbsnp
          type: File
        - name: genome_resources__variation__cosmic
          type: File
        - name: reference__genome_context
          type:
          - 'null'
          - string
          - items:
            - 'null'
            - string
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
        - name: genome_resources__aliases__ensembl
          type: string
        - name: reference__rtg
          type: File
        - name: genome_resources__aliases__snpeff
          type: string
        - name: align_bam
          type:
          - File
          - 'null'
        - name: regions__sample_callable
          type:
          - File
          - 'null'
        - name: config__algorithm__callable_regions
          type: File
        - name: vrn_file_joint
          type: File
        name: jointvc_rec
        type: record
      type: array
    type: array
baseCommand:
- bcbio_nextgen.py
- runfn
- summarize_vc
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=variants__calls,variants__gvcf,variants__samples,validate__grading_summary,validate__grading_plots
- sentinel_inputs=jointvc_rec:record
outputs:
- id: variants__calls
  type:
    items:
      items:
      - File
      - 'null'
      type: array
    type: array
- id: variants__gvcf
  type:
    items:
    - 'null'
    - items:
      - File
      - 'null'
      type: array
    type: array
- id: variants__samples
  type:
    items:
      items:
        items:
        - File
        - 'null'
        type: array
      type: array
    type: array
- id: validate__grading_summary
  type:
    items:
    - File
    - 'null'
    type: array
- id: validate__grading_plots
  type:
    items:
      items:
      - File
      - 'null'
      type: array
    type: array
$namespaces:
  dx: https://www.dnanexus.com/cwl#
