cwlVersion: v1.2
class: Workflow

inputs:
  file1:
    type: File
    default:
      class: File
      location: dx://file-G0G0V1Q0yzZZZXfx3xPK1B1Z

outputs:
  o:
    type: File
    outputSource: step1/o

steps:
  step1:
    in:
      catfile1: file1
    out: [o]
    run: 
      class: CommandLineTool

      inputs:
        catfile1:
          type: File

      outputs:
        o:
          type: File
          outputBinding: { glob: output }

      arguments: [cat, $(inputs.catfile1.path)]
      stdout: output
