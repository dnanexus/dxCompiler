#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
stdout: output.txt
inputs:
  example_in:
    type:
    - 'null'
    - type: enum
      symbols:
      - string1
      - string2
    inputBinding:
      prefix: --example-string
      position: 0
baseCommand: [echo]
outputs:
  example_out:
    type: stdout