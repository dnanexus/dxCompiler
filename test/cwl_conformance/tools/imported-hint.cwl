#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
inputs: []
outputs:
  out: stdout

hints:
- $import: envvar.yml

baseCommand: ["/bin/sh", "-c", "echo $TEST_ENV"]

stdout: out
