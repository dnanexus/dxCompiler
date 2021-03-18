#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

inputs:
  file1:
    type: File

outputs:
  count_output:
    type: int
    outputSource: step2/output

steps:
  step1:
    run:
      class: CommandLineTool
      inputs:
        file1: File
      outputs:
        output:
          type: File
          outputBinding: { glob: output }
      baseCommand: [ wc, -l ]
      stdin: $(inputs.file1.path)
      stdout: output
    in:
      file1: file1
    out: [output]

  step2:
    run:
      class: ExpressionTool
      requirements:
        - class: InlineJavascriptRequirement
      cwlVersion: v1.2
      inputs:
        file1:
          type: File
          loadContents: true
      outputs:
        output: int
      expression: "$({'output': parseInt(inputs.file1.contents)})"
    in:
      file1: step1/output
    out: [output]
