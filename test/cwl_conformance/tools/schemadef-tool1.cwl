#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

requirements:
  - class: SchemaDefRequirement
    types:
      - name: HelloType
        type: record
        fields:
          - name: a
            type: string
          - name: b
            type: string

inputs:
    - id: hello
      type: HelloType
      inputBinding:
        valueFrom: $(self.a)/$(self.b)

outputs:
    - id: output
      type: File
      outputBinding:
        glob: output.txt

stdout: output.txt
baseCommand: echo
