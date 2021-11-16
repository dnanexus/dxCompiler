#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  class: Workflow
  inputs:
  - id: jamie
    type: File
  outputs:
  - id: md5
    outputSource: "sum/md5"
    type: string
  steps:
  - id: sum
    in:
    - id: jamie
      source: "jamie"
    out:
    - id: md5
    requirements:
      DockerRequirement:
        dockerPull: "ubuntu:latest"
    run:
      inputs:
      - id: jamie
        type: File
      outputs:
      - id: md5
        outputBinding:
          glob: md5-stdOut.txt
          loadContents: true
          outputEval: $(self[0].contents)
        type: string
      class: CommandLineTool
      requirements:
      - class: ShellCommandRequirement
      - class: InlineJavascriptRequirement
      arguments:
      - valueFrom: "/usr/bin/md5sum"
        shellQuote: false
      - valueFrom: $(inputs.jamie)
        shellQuote: true
      - valueFrom: '|'
        shellQuote: false
      - valueFrom: cut
        shellQuote: false
      - valueFrom: -c1-32
        shellQuote: false
      stdout: md5-stdOut.txt
      hints:
        NetworkAccess:
          networkAccess: true
        LoadListingRequirement:
          loadListing: deep_listing
