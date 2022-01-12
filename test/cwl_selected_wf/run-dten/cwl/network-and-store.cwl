#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: network-and-store
label: network-and-store
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
steps:
  run-network:
    in:
      beta: beta
      mu: mu
      w: w
      protein-list: protein-list
      condition: condition
    out: network-file
    run: steps/run-network-with-params.cwl
  store-network:
    in:
      file_to_store: run-network/network-file
      parentid: output-folder-id
      synapse_config: synapse_config
    run: cwl/synapse-store-tool.cwl
    out: []
outputs:
- id: network-file
  type: File
  outputSource: run-network/network-file
