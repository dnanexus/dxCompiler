#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: store-tables
label: store-tables
requirements:
- class: DockerRequirement
  dockerPull: sage-bionetworks/synapse
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  nodeTable:
    type: File
  termTable:
    type: File
  output-project-id:
    type: string
  synapse_config:
    type: File
baseCommand: ['synapse store']
outputs: []
