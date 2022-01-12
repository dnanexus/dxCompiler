#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: run-meta-analysis
label: run-meta-analysis
requirements:
- class: DockerRequirement
  dockerPull: sgosline/dten
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(inputs.synapse_config)
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  synapse_config:
    type: File
  input:
    type: File[]
    inputBinding:
      prefix: '-i'
      itemSeparator: ","
  output:
    type: string
    inputBinding:
      prefix: '-o'
  project:
    type: string
    inputBinding:
      prefix: '-p'
  folder:
    type: string
    inputBinding:
      prefix: '-f'
baseCommand: ['Rscript', '/usr/local/bin/metaNetworkComparisons.R']
outputs:
  nodefile:
    type: File
    outputBinding:
      glob: "*nodeOutput.tsv"
  termfile:
    type: File
    outputBinding:
      glob: "*termOutput.tsv"
