#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
id: get-mv-samples
label: get-mv-samples
requirements:
- class: DockerRequirement
  dockerPull: sgosline/dten
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(inputs.synapse_config)
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  synapse_config:
    type: File
  parent-folder:
    type: string
    inputBinding:
      prefix: "--folderid"
  fileview:
    type: string
    inputBinding:
      prefix: "--fileview"
  num-samps:
    type: string
    inputBinding:
      prefix: "--n"
baseCommand: ['Rscript', '/usr/local/bin/subSampleTumorTypes.R']
stdout: message
outputs:
  synids:
    type: string[]
    outputBinding:
      glob: message
      loadContents: true
      outputEval: $(self[0].contents.split('\n'))
