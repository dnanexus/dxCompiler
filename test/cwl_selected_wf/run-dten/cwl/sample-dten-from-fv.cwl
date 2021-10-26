label: sample-dten-from-fv
id: sample-dten-from-fv
cwlVersion: v1.0
class: Workflow

inputs:
  input-fileview:
    type: string
  metaviper-folder:
    type: string
  synapse-config:
    type: File
  num-samps:
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
  get-prots-from-query:
    run: steps/get-mv-samples.cwl
    in:
      synapse_config: synapse-config
      parent-folder: metaviper-folder
      fileview: input-fileview
      num-samps: num-samps
    out:
      [synids]
  get-prot-files:
    run: https://raw.githubusercontent.com/Sage-Bionetworks-Workflows/cwl-tool-synapseclient/main/cwl/synapse-get-tool.cwl
    scatter: synapseid
    in:
      synapseid: get-prots-from-query/synids
      synapse_config: synapse-config
    out:
      [filepath]
  build-networks:
    scatter: [beta, mu, w]
    scatterMethod: flat_crossproduct
    run: build-store-networks-with-params.cwl
    in:
      beta: beta-params
      mu: mu-params
      w: w-params
      protein-lists: get-prot-files/filepath
      output-project-id: output-project-id
      output-folder-id: output-parent-id
      synapse_config: synapse-config
      net-name: name
    out:
      [network-file]
