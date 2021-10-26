#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: proteins-from-genes
label: proteins-from-genes
requirements:
- class: DockerRequirement
  dockerPull: sgosline/dten
- class: InlineJavascriptRequirement
hints:
  NetworkAcess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  gene-data:
    type: File[]
    inputBinding:
      prefix: '-i'
      itemSeparator: ','
  id-type:
    type: string
    inputBinding:
      prefix: '-d'
  condition:
    type: string[]?
    inputBinding:
      prefix: '-c'
      itemSeparator: ','
baseCommand: ['Rscript', '/usr/local/bin/runMetaViper.R']
outputs:
  protein-lists:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.tsv"
