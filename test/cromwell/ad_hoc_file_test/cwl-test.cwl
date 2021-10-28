#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: example.sh
      entry: |-
        PREFIX='Message is:'
        MSG="\${PREFIX} Hello world!"
        echo \${MSG}
hints:
  DockerRequirement:
    dockerPull: ubuntu:latest
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs: []
baseCommand: ["sh", "example.sh"]
stdout: output.txt
outputs:
  example_out:
    type: stdout
