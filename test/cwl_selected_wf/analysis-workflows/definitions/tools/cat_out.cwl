#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:xenial"
- class: ResourceRequirement
  ramMin: 4000
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  pindel_outs:
    type: File[]
    inputBinding:
      position: 1
baseCommand: ['/bin/cat']
stdout: "per_chromosome_pindel.out"
outputs:
  pindel_out:
    type: stdout
