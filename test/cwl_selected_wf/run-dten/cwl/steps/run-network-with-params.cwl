label: run-network-with-params
id: run-network-with-params
cwlVersion: v1.0
class: CommandLineTool

baseCommand: ['Rscript','/usr/local/bin/runNetworkFromGenes.R']


requirements:
  - class: DockerRequirement
    dockerPull: sgosline/dten
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.synapse_config)

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

outputs:
  network-file:
    type: File
    outputBinding:
      glob: "*.rds"
