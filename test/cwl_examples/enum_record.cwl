#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo

inputs:
  message:
    type: string
    inputBinding:
      position: 1
  second:
    type:
      - type: enum
        symbols: [homo_sapiens, mus_musculus]
  first:
    type:
      type: record
      fields:
        species:
          - type: enum
            symbols: [homo_sapiens, mus_musculus]
  file1:
    type: File
    label: Input File
    doc: "The file that will be copied using 'cat'"

  args.py:
    type: File
    default:
      class: File
      location: args.py
    inputBinding:
      position: -1
outputs:
  - id: sam
    type: [File]
    outputBinding: { glob: output.sam }