#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
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
  region_pindel_outs:
    type: File[]
    inputBinding:
      position: -1
baseCommand: ['/bin/cat']
arguments: [{shellQuote: false, valueFrom: "|"}, "/bin/grep", "ChrID", "/dev/stdin"]
stdout: "all_region_pindel.head"
outputs:
  all_region_pindel_head:
    type: stdout
