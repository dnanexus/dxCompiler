#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool
doc: ["hello", "world"]
inputs:
  - id: reference
    type: string
    inputBinding: { position: 1 }
outputs:
  - id: sam
    type: File
    outputBinding: { glob: output.sam }
