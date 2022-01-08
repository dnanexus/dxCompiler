#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: MPLCONFIGDIR
    envValue: .
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
hints: []
inputs:
- id: batch_rec
  type:
    items:
      fields:
      - name: description
        type: string
      - name: resources
        type: string
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
      - name: vrn_file
        type:
        - 'null'
        - string
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
      name: batch_rec
      type: record
    type: array
steps:
- id: get_parallel_regions
  in:
  - id: batch_rec
    source: batch_rec
  out:
  - id: region_block
  run: steps/get_parallel_regions.cwl
- id: variantcall_batch_region
  in:
  - id: batch_rec
    source: batch_rec
  - id: region_block_toolinput
    source: get_parallel_regions/region_block
  out:
  - id: vrn_file_region
  - id: region_block
  run: steps/variantcall_batch_region.cwl
  scatter:
  - region_block_toolinput
  scatterMethod: dotproduct
- id: concat_batch_variantcalls
  in:
  - id: batch_rec
    source: batch_rec
  - id: region_block
    source: variantcall_batch_region/region_block
  - id: vrn_file_region
    source: variantcall_batch_region/vrn_file_region
  out:
  - id: vrn_file
  run: steps/concat_batch_variantcalls.cwl
- id: compare_to_rm
  in:
  - id: batch_rec
    source: batch_rec
  - id: vrn_file
    source: concat_batch_variantcalls/vrn_file
  out:
  - id: vc_rec
  run: steps/compare_to_rm.cwl
outputs:
- id: vc_rec
  outputSource: compare_to_rm/vc_rec
  type:
    items:
      fields:
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
      - name: description
        type: string
      - name: resources
        type: string
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
      name: vc_rec
      type: record
    type: array
