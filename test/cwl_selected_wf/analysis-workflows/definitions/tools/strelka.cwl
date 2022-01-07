#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "strelka 2.9.9"
requirements:
  ResourceRequirement:
    coresMin: 4
    ramMin: 4000
  DockerRequirement:
    dockerPull: "mgibio/strelka-cwl:2.9.9"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  tumor_bam:
    type: File
    inputBinding:
      prefix: '--tumorBam='
      separate: false
      position: 3
    secondaryFiles: [.bai, ^.bai]
  normal_bam:
    type: File
    inputBinding:
      prefix: '--normalBam='
      separate: false
      position: 4
    secondaryFiles: [.bai, ^.bai]
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: '--referenceFasta='
      separate: false
      position: 5
  exome_mode:
    type: boolean
    inputBinding:
      prefix: '--exome'
      position: 6
  cpu_reserved:
    type: int?
baseCommand: ["/usr/bin/perl", "/usr/bin/docker_helper.pl"]
arguments: [{valueFrom: $(inputs.cpu_reserved), position: 1}, {valueFrom: $(runtime.outdir),
    position: 2}]
outputs:
  indels:
    type: File
    outputBinding:
      glob: "results/variants/somatic.indels.vcf.gz"
  snvs:
    type: File
    outputBinding:
      glob: "results/variants/somatic.snvs.vcf.gz"
