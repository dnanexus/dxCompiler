#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: "synapse-get"
label: "Synapse Get Tool"
requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: .synapseConfig
      entry: $(inputs.synapse_config)
hints:
  DockerRequirement:
    dockerPull: sagebionetworks/synapsepythonclient:v2.4.0
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
- id: synapse_config
  type: File
- id: synapseid
  type: string
baseCommand: synapse
arguments:
- valueFrom: get
- valueFrom: $(inputs.synapseid)
outputs:
- id: filepath
  type: File
  outputBinding:
    glob: '*'
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-5841-0198
  s:email: thomas.yu@sagebionetworks.org
  s:name: Thomas Yu
