#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: "synapse-set-annotations"
label: "Synapse set annotations tool"
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: .synapseConfig
    entry: $(inputs.synapse_config)
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
- id: synapseid
  type: string
  inputBinding:
    prefix: --id
- id: annotations_text
  type: string
  inputBinding:
    prefix: --annotations
- id: replace
  type: boolean
  default: true
  inputBinding:
    prefix: --replace
baseCommand:
- synapse
- set-annotations

outputs: []
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-0326-7494
  s:email: andrew.lamb@sagebase.org
  s:name: Andrew Lamb
