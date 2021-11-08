#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: cwl_stdout_expression
  class: CommandLineTool
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
  inputs:
  - id: foo
    type: string
  - id: bar
    type: string
  outputs:
    b: stdout
  stdout: "stdout-$(inputs.foo)-$(inputs.bar).txt"
  baseCommand: []
  arguments: [echo, $(inputs.foo), $(inputs.bar)]
