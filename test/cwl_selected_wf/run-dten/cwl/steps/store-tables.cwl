label: store-tables
id: store-tables
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: sage-bionetworks/synapse

baseCommand: ['synapse store']

inputs:
  nodeTable:
    type: File
  termTable:
    type: File
  output-project-id:
    type: string
  synapse_config:
    type: File

outputs:
  []
