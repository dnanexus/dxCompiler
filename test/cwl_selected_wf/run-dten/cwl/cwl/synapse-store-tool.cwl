#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: "synapse-store"
label: "Synapse command line client subcommand for storing a file."
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
- id: file_to_store
  type: File
- id: parentid
  type: string
- id: name
  type: string?
- id: used
  type: string[]?
- id: executed
  type: string[]?
baseCommand: synapse
arguments:
- valueFrom: store
- valueFrom: $(inputs.parentid)
  prefix: --parentId
- valueFrom: $(inputs.used)
  prefix: --used
- valueFrom: $(inputs.executed)
  prefix: --executed
- valueFrom: $(inputs.name)
  prefix: --name
- valueFrom: $(inputs.file_to_store.path)
  prefix: --
stdout: stdout.txt
outputs:
- id: stdout
  type: File
  outputBinding:
    glob: stdout.txt
- id: file_id
  type: string
  outputBinding:
    glob: stdout.txt
    loadContents: true
    outputEval: $(self[0].contents.split("\n")[5].split(/(\s+)/)[4])
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0001-5729-7376
  s:email: kenneth.daily@sagebionetworks.org
  s:name: Kenneth Daily
s:contributor:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-0326-7494
  s:email: andrew.lamb@sagebase.org
  s:name: Andrew Lamb
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-5841-0198
  s:email: thomas.yu@sagebase.org
  s:name: Thomas Yu
