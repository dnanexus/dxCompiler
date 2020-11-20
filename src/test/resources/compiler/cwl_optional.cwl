#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo

inputs:
  first:
    type:
      - "null"
      - type: enum
        symbols: [homo_sapiens, mus_musculus]

outputs:
  - id: sam
    type: ["null",  File]
    outputBinding: { glob: output.sam }