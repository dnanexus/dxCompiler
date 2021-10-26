#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
id: run-dten
label: run-dten
requirements:
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: MultipleInputFeatureRequirement
inputs:
  conditions:
    type: string[]
  queries:
    type: string[]
  synapse-config:
    type: File
  gene-id-type:
    type: string
  beta-params:
    type: double[]
  mu-params:
    type: double[]
  w-params:
    type: double[]
  output-parent-id:
    type: string
  output-project-id:
    type: string
  metaviper-store-id:
    type: string
  name:
    type: string
steps:
  download-file:
    run: cwl/synapse-query-tool.cwl
    scatter: query
    in:
      query: queries
      synapse_config: synapse-config
    out: [query_result]
  get-prots:
    run: steps/proteins-from-genes.cwl
    in:
      gene-data: download-file/query_result
      condition: conditions
      id-type: gene-id-type
    out: [protein-lists]
  store-prots:
    run: cwl/synapse-store-tool.cwl
    scatter: file_to_store
    in:
      file_to_store: get-prots/protein-lists
      synapse_config: synapse-config
      parentid: metaviper-store-id
    out: []
  build-networks:
    scatter: [beta, mu, w]
    scatterMethod: flat_crossproduct
    run: build-store-networks-with-params.cwl
    in:
      beta: beta-params
      mu: mu-params
      w: w-params
      protein-lists: get-prots/protein-lists
      output-project-id: output-project-id
      output-folder-id: output-parent-id
      synapse_config: synapse-config
      net-name: name
    out: [network-file]
outputs:
  out:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: build-networks/network-file
