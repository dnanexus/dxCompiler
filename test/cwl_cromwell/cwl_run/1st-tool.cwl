#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  message:
    type: string
    inputBinding:
      position: 1
baseCommand: echo
outputs: []
