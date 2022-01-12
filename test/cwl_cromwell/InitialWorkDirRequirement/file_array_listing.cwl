#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  class: CommandLineTool
  baseCommand: ['ls']
  requirements:
  - class: DockerRequirement
    dockerPull: "python:3.5.0"
  - class: InitialWorkDirRequirement
    listing: $(inputs.files)
  - class: InlineJavascriptRequirement
  stdout: file_list
  inputs:
  - id: files
    type:
      type: array
      items: File
    default: ["dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:/test_data/cwl/cwl_secondary_files/foo.txt"]
  arguments: ["*.txt"]
  outputs:
  - id: file_list
    type: string
    outputBinding:
      glob: file_list
      loadContents: true
      outputEval: $(self[0].contents.trim())
  hints:
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
