#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
inputs: []
outputs:
  out: stdout

hints:
- class: EnvVarRequirement
  envDef:
    - envName: "TEST_ENV"
      envValue: "hello test env"

baseCommand: ["/bin/sh", "-c", "echo $TEST_ENV"]

stdout: out
