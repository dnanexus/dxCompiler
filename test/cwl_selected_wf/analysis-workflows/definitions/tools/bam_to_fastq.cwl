#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Picard: BAM to FASTQ"
requirements:
- class: ResourceRequirement
  coresMin: 1
  ramMin: 6000
  tmpdirMin: 25000
- class: DockerRequirement
  dockerPull: "mgibio/rnaseq:1.0.0"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  bam:
    type: File
    inputBinding:
      prefix: "I="
      separate: false
      position: 1
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/picard/picard.jar", "SamToFastq",
  "VALIDATION_STRINGENCY=SILENT"]
arguments: [valueFrom: "F=$(runtime.outdir)/read1.fastq", valueFrom: "F2=$(runtime.outdir)/read2.fastq"]
outputs:
  fastq1:
    type: File
    outputBinding:
      glob: "read1.fastq"
  fastq2:
    type: File
    outputBinding:
      glob: "read2.fastq"
