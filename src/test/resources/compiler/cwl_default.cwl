#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo

inputs:
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