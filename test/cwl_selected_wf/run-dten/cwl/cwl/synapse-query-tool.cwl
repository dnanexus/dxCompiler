#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: synapse-query-tool
label: synapse-query-tool
requirements:
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
  synapse_config:
    type: File
    doc: synapseConfig file
  query:
    type: string
    inputBinding:
      position: 1
      prefix: query
    doc: query
baseCommand: synapse
stdout: query_result.tsv
outputs:
- id: query_result
  type: stdout
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-5841-0198
  s:email: thomas.yu@sagebionetworks.org
  s:name: Thomas Yu
