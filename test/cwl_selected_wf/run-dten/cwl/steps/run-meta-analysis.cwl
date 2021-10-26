label: run-network-with-params
id: run-network-with-params
cwlVersion: v1.0
class: CommandLineTool

baseCommand: ['Rscript','/usr/local/bin/metaNetworkComparisons.R']


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

outputs:
  nodefile:
    type: File
    outputBinding:
      glob: "*nodeOutput.tsv"
  termfile:
    type: File
    outputBinding:
      glob: "*termOutput.tsv"
