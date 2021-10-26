#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: "synapse-get-sts"
label: "Synapse Get STS Token Tool"
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
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
- id: permission
  type: string
baseCommand: synapse
arguments:
- valueFrom: get-sts-token
- valueFrom: $(inputs.synapseid)
- valueFrom: $(inputs.permission)
- valueFrom: json
  prefix: --output
stdout: output.json
outputs:
- id: json_out
  type: stdout
- id: bucket
  type: string
  outputBinding:
    glob: output.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['bucket'])
- id: basekey
  type: string
  outputBinding:
    glob: output.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['baseKey'])
- id: accesskey_id
  type: string
  outputBinding:
    glob: output.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['accessKeyId'])
- id: secret_accesskey
  type: string
  outputBinding:
    glob: output.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['secretAccessKey'])
- id: session_token
  type: string
  outputBinding:
    glob: output.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['sessionToken'])
$namespaces:
  s: https://schema.org/
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-5841-0198
  s:email: thomas.yu@sagebionetworks.org
  s:name: Thomas Yu
