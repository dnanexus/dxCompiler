#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo

inputs:
  second:
    type:
      - type: enum
        symbols: [homo_sapiens, mus_musculus]
outputs:
  - id: sam
    type: [File]
    outputBinding: { glob: output.sam }