label: run-dten
id: run-dten
cwlVersion: v1.0
class: Workflow

inputs:
  input-query-list:
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

outputs:
  out:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: build-networks/network-file

requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement

steps:
  download-file:
    run: https://raw.githubusercontent.com/Sage-Bionetworks-Workflows/cwl-tool-synapseclient/main/cwl/synapse-query-tool.cwl
    scatter: query
    in:
      query: input-query-list
      synapse_config: synapse-config
    out: [query_result]
  get-prots:
    run: steps/proteins-from-genes.cwl
    in:
      gene-data: download-file/query_result
      id-type: gene-id-type
    out:
      [protein-lists]
  store-prots:
    run: https://raw.githubusercontent.com/Sage-Bionetworks-Workflows/cwl-tool-synapseclient/main/cwl/synapse-store-tool.cwl
    scatter: file_to_store
    in:
      file_to_store:  get-prots/protein-lists
      synapse_config: synapse-config
      parentid: metaviper-store-id
    out:
      []
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
    out:
      [network-file]
