#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: "ubuntu:latest"
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  message:
    type: string[]?
  unsupplied_optional:
    type: string[]?
baseCommand: echo
outputs: []
