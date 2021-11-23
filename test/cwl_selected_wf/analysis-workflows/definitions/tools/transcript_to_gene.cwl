#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Kallisto: TranscriptToGene"
requirements:
- class: ResourceRequirement
  ramMin: 2000
  coresMin: 1
- class: DockerRequirement
  dockerPull: "mgibio/rnaseq:1.0.0"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  gene_transcript_lookup_table:
    type: File
    inputBinding:
      position: 1
  transcript_table_h5:
    type: File
    inputBinding:
      position: 2
baseCommand: ["/usr/local/bin/Rscript"]
arguments: ["/usr/src/transcript_to_gene.R"]
outputs:
  gene_abundance:
    type: File
    outputBinding:
      glob: "gene_abundance.tsv"
