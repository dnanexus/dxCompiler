#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement:
    expressionLib:
    - "function foo() { return 2; }"
hints:
  DockerRequirement:
    dockerPull: "ubuntu:latest"
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs: []
arguments: [echo, $(foo())]
stdout: whatever.txt
outputs:
  out: stdout
