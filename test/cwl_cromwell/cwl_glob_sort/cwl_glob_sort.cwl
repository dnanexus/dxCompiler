#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  cwlVersion: v1.0
  class: CommandLineTool
  requirements:
  - class: InlineJavascriptRequirement
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
  inputs: []
  baseCommand: [touch, z, y, x, w, c, b, a]
  outputs:
    letters:
      type: string
      outputBinding:
        glob: '?'
        outputEval: |
          ${ return self.sort(function(a,b) { return a.location > b.location ? 1 : (a.location < b.location ? -1 : 0) }).map(f => f.basename).join(" ") }
