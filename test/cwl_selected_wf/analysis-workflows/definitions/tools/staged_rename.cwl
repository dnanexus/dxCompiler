#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Staged Renamer"
doc: "Renames a file by staging and then `mv`ing it.  A workaround for workflow engines\
  \ that don't support rename.cwl.  If running in cwltool, use the other one instead."
requirements:
- class: ResourceRequirement
  ramMin: 4000
  coresMin: 1
- class: DockerRequirement
  dockerPull: ubuntu:bionic
- class: InitialWorkDirRequirement
  listing:
  - $(inputs.original)
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  original:
    type: File
  name:
    type: string
baseCommand: ["/bin/mv"]
arguments: [$(inputs.original.basename), $(inputs.name)]
outputs:
  replacement:
    type: File
    outputBinding:
      glob: $(inputs.name)
