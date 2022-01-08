#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  class: CommandLineTool
  cwlVersion: v1.0
  doc: "Asks for CPU and Memory minimums"
  requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 2
    coresMax: 2
    ramMin: 6675.72
    ramMax: 6675.72
  hints:
    DockerRequirement:
      dockerPull: python:latest
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
  inputs: []
  outputs:
    machine_type:
      type: string
      outputBinding:
        glob: stdout
        loadContents: true
        outputEval: $(self[0].contents.trim())
  baseCommand: ['curl', 'http://metadata.google.internal/computeMetadata/v1/instance/machine-type',
    '-H', 'Metadata-Flavor: Google']
  stdout: stdout
