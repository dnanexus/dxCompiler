#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
inputs: []
steps:
  step1:
    run: tools/task.cwl
    in:
      example_in: {default: string1}
    out: [example_out]
outputs:
  wf_out:
    type: File
    outputSource: step1/example_out
