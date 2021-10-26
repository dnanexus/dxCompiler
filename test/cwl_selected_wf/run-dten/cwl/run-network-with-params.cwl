#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
id: run-network-with-params
label: run-network-with-params
requirements:
- class: DockerRequirement
  dockerPull: sgosline/dten
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(inputs.synapse_config)
hints:
  NetworkAcess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  synapse_config:
    type: File
  mu:
    type: double
    inputBinding:
      prefix: '-m'
  beta:
    type: double
    inputBinding:
      prefix: '-b'
  w:
    type: double
    inputBinding:
      prefix: '-w'
  protein-list:
    type: File
    inputBinding:
      prefix: '-i'
baseCommand: ['Rscript', '/usr/local/bin/runNetworkFromGenes.R']
outputs:
  network-file:
    type: File
    outputBinding:
      glob: "*.rds"
