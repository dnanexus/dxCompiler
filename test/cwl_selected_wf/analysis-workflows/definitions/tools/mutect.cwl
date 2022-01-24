#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Mutect2 (GATK 4)"
requirements:
- class: ResourceRequirement
  ramMin: 32000
  tmpdirMin: 100000
- class: DockerRequirement
  dockerPull: "broadinstitute/gatk:4.1.8.1"
- class: InitialWorkDirRequirement
  listing:
  - entryname: 'Mutect2.sh'
    entry: |
      set -o pipefail
      set -o errexit

      export tumor_bam="$3"
      export normal_bam="$4"

      NORMAL=`samtools view -H $normal_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
      TUMOR=`samtools view -H $tumor_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`

      /gatk/gatk Mutect2 --java-options "-Xmx20g" -O $1 -R $2 -I $3 -tumor "$TUMOR" -I $4 -normal "$NORMAL" -L $5 #Running Mutect2.
      /gatk/gatk FilterMutectCalls -R $2 -V mutect.vcf.gz -O mutect.filtered.vcf.gz #Running FilterMutectCalls on the output vcf.
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      position: 2
  tumor_bam:
    type: File
    inputBinding:
      position: 3
    secondaryFiles: [.bai]
  normal_bam:
    type: File?
    inputBinding:
      position: 4
    secondaryFiles: [.bai]
  interval_list:
    type: File
    inputBinding:
      position: 5
baseCommand: ["/bin/bash", "Mutect2.sh"]
arguments:
- position: 1
  valueFrom: mutect.vcf.gz
outputs:
  vcf:
    type: File
    outputBinding:
      glob: "mutect.filtered.vcf.gz"
    secondaryFiles: [.tbi]
