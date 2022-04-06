#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
requirements:
- class: ScatterFeatureRequirement
- class: DockerRequirement
  dockerPull: ubuntu:latest
inputs:
  inp: string[]
  inp2: string
steps:
  step1:
    in:
      echo_in: inp
      echo_in2: inp2
    out: [echo_out]
    scatter: echo_in
    run:
      class: CommandLineTool
      inputs:
        echo_in:
          type: string
          inputBinding: {}
        echo_in2:
          type: string
          inputBinding: {}
      outputs:
        echo_out:
          type: string
          outputBinding:
            glob: "step1_out"
            loadContents: true
            outputEval: $(self[0].contents)
      baseCommand: "echo"
      arguments:
      - "-n"
      - "foo"
      stdout: "step1_out"
      hints:
        NetworkAccess:
          networkAccess: true
        LoadListingRequirement:
          loadListing: deep_listing
outputs:
  out:
    type: string[]
    outputSource: step1/echo_out
