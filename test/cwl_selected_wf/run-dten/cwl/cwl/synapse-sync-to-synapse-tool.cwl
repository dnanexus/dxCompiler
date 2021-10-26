#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: "sync-to-synapse"
label: "Sync to Synapse tool"
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: .synapseConfig
    entry: $(inputs.synapse_config)
  - $(inputs.files)

hints:
  DockerRequirement:
    dockerPull: sagebionetworks/synapsepythonclient:v2.4.0
  NetworkAcess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
- id: synapse_config
  type: File
- id: files
  type: File[]
- id: manifest_file
  type: File
baseCommand: synapse
arguments:
- valueFrom: sync
- valueFrom: $(inputs.manifest_file.path)
outputs: []
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-0326-7494
  s:email: andrew.lamb@sagebase.org
  s:name: Andrew Lamb
s:contributor:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-5841-0198
  s:email: thomas.yu@sagebionetworks.org
  s:name: Thomas Yu
