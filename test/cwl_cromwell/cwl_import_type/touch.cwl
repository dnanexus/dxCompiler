#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: DockerRequirement
  dockerPull: ubuntu:bionic-20180426
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
- id: input
  type: string
  inputBinding:
    position: 0
baseCommand: [touch]
outputs:
- id: output
  type: File
  outputBinding:
    glob: $(inputs.input)
