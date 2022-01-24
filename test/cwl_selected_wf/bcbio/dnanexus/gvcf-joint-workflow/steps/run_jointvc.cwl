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
  ramMin: 2560
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: gatk4
    specs:
    - https://anaconda.org/bioconda/gatk4
  - package: gatk
    specs:
    - https://anaconda.org/bioconda/gatk
  - package: gvcftools
    specs:
    - https://anaconda.org/bioconda/gvcftools
  - package: sentieon
    specs:
    - https://anaconda.org/bioconda/sentieon
- class: arv:APIRequirement
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing
inputs:
- id: jointvc_batch_rec
  type:
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
      name: jointvc_batch_rec
      type: record
    type: array
- id: region_toolinput
  type: string
baseCommand:
- bcbio_nextgen.py
- runfn
- run_jointvc
- cwl
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=batch-parallel
- sentinel_outputs=vrn_file_region,region
- sentinel_inputs=jointvc_batch_rec:record,region:var
outputs:
- id: vrn_file_region
  secondaryFiles:
  - .tbi
  type:
  - File
  - 'null'
- id: region
  type: string
$namespaces:
  arv: http://arvados.org/cwl#
  dx: https://www.dnanexus.com/cwl#
