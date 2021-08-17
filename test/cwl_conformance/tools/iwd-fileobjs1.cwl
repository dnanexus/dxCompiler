#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
doc: Can have a file declared directly in InitialWorkDir
requirements:
  InitialWorkDirRequirement:
    listing:
      - class: File
        location: "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:/test_data/cwl/loadContents/inp-filelist.txt"
      - class: Directory
        location: "dx://project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq:/test_data/cwl/testdir/"
inputs: []
outputs:
  filelist:
    type: File
    outputBinding:
      glob: inp-filelist.txt
  testdir:
    type: Directory
    outputBinding:
      glob: testdir
baseCommand: "true"
