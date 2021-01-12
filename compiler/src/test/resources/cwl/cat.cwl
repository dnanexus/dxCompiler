cwlVersion: v1.2
class: CommandLineTool
doc: Write a file to stdout using cat
baseCommand: cat
stdout: output.txt
inputs:
  file:
    type: File
    inputBinding:
      position: 1
outputs:
  contents:
      type: stdout