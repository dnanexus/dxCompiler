#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: "synapse-create"
label: "Synapse command line client subcommand for creating a project/folder."
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
- id: parentid
  type: string?
  inputBinding:
    prefix: --parentId
- id: name
  type: string
  inputBinding:
    prefix: --name
- id: description
  type: string?
  inputBinding:
    prefix: --description
- id: description_file
  type: File?
  inputBinding:
    prefix: --descriptionFile
- id: type
  type:
    type: enum
    symbols: [Project, Folder]
  inputBinding:
    position: 2
baseCommand:
- synapse
- create

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
    outputEval: $(self[0].contents.split("\n")[0].split(/(\s+)/)[4])
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-4621-1589
  s:email: bruno.grande@sagebase.org
  s:name: Bruno Grande
