#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

inputs:
  special_file:
    type: File
    default:
      class: File
      location: dx://file-G1353z00yzZYkb2V4j4g7yb9

outputs:
  cores:
    type: File
    outputSource: report/output

steps:
  count:
    in:
      special_file: special_file
    out: [output]
    run: dynresreq.cwl

  report:
    in:
      file1: count/output
    out: [output]
    run: cat-tool.cwl
