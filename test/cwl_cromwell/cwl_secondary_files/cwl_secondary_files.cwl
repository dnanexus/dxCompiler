#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  class: CommandLineTool
  requirements:
    - class: InlineJavascriptRequirement
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
  inputs:
  - id: command
    type: string
  - id: f
    type: File
    inputBinding:
      position: 2
    secondaryFiles: [.also]
  - id: of
    type: File?
    inputBinding:
      position: 3
    secondaryFiles: [.also]
  - id: fs
    type:
      type: array
      items: File
      inputBinding:
        position: 4
    secondaryFiles: [.also]
  - id: fr
    secondaryFiles: [.also]
    type:
      type: record
      name: foo
      fields:
      - name: a
        type: File
        inputBinding:
          position: 5
  outputs:
    the_answer:
      type: string
      outputBinding:
        outputEval: ${ return "$(" + 42 + ")"; }
  baseCommand: []
  arguments: ["bash", "-c", $(inputs.command)]
