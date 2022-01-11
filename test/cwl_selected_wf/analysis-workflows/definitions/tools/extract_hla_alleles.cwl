#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "ubuntu:xenial"
- class: ResourceRequirement
  ramMin: 2000
- class: ShellCommandRequirement
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  allele_file:
    type: File
    inputBinding:
      position: -1
baseCommand: ['/usr/bin/awk', '{getline; printf "HLA-"$2 " HLA-"$3 " HLA-"$4 " HLA-"$5
    " HLA-"$6 " HLA-"$7}']
arguments: [{shellQuote: false, valueFrom: ">"}, 'helper.txt']
outputs:
  allele_string:
    type: string[]
    outputBinding:
      glob: helper.txt
      loadContents: true
      outputEval: $(self[0].contents.split(" "))
