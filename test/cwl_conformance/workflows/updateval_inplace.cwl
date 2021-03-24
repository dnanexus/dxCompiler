#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.r)
        writable: true
  InplaceUpdateRequirement:
    inplaceUpdate: true
hints:
  DockerRequirement:
    dockerPull: "python:2.7.15-alpine3.7"
inputs:
  r: File
  script:
    type: File
    default:
      class: File
      location: dx://file-G0G3BZj0yzZg3BBz3x4Y2Qkq
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.r.basename)
arguments: [python, $(inputs.script), $(inputs.r.basename)]
