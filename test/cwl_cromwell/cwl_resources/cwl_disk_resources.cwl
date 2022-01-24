#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$namespaces:
  dx: https://www.dnanexus.com/cwl#
$graph:
- id: main
  class: CommandLineTool
  cwlVersion: v1.0
  doc: "Asks for disk minimums"
  requirements:
    - class: ResourceRequirement
      outdirMin: 10240
      tmpdirMin: 10240
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
  hints:
  - class: DockerRequirement
    dockerPull: gcr.io/google.com/cloudsdktool/cloud-sdk:slim
  - class: dx:InputResourceRequirement
    indirMin: 10240
  - class: NetworkAccess
    networkAccess: true
  - class: LoadListingRequirement
    loadListing: deep_listing
  inputs: []
  outputs:
    disk_size:
      type: int
      outputBinding:
        glob: stdout
        loadContents: true
        outputEval: $(parseInt(self[0].contents.trim()))
  arguments:
  - valueFrom: |
      apt-get install --assume-yes jq > /dev/null && NAME=`curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/name` && ZONE=`basename \`curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone\`` && PROJECT=`curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/project/project-id` && curl -H "Authorization: Bearer `gcloud auth print-access-token`" "https://www.googleapis.com/compute/v1/projects/${PROJECT}/zones/${ZONE}/disks/${NAME}-1?fields=sizeGb" | jq -r '.sizeGb'
    shellQuote: false
  stdout: stdout
