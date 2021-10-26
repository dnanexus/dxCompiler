label: network-and-store
id: network-and-store
cwlVersion: v1.0
class: Workflow


inputs:
  beta:
    type: double
  mu:
    type: double
  w:
    type: double
  protein-list:
    type: File
  condition:
    type: string
  output-folder-id:
    type: string
  synapse_config:
    type: File

outputs:
  - id: network-file
    type: File
    outputSource: run-network/network-file

steps:
  run-network:
    in:
      beta: beta
      mu: mu
      w: w
      protein-list: protein-list
      condition: condition
    out:
      network-file
    run: steps/run-network-with-params.cwl
  store-network:
    in:
      file_to_store: run-network/network-file
      parentid: output-folder-id
      synapse_config: synapse_config
    run: https://raw.githubusercontent.com/Sage-Bionetworks-Workflows/cwl-tool-synapseclient/main/cwl/synapse-store-tool.cwl
    out:
      []
