#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "expand interval list regions by a given number of basepairs"
requirements:
- class: ResourceRequirement
  ramMin: 4000
- class: DockerRequirement
  dockerPull: "broadinstitute/picard:2.23.6"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  interval_list:
    type: File
    inputBinding:
      prefix: "INPUT="
      separate: false
  roi_padding:
    type: int
    inputBinding:
      prefix: "PADDING="
      separate: false
baseCommand: ["/usr/bin/java", "-Xmx3g", "-jar", "/usr/picard/picard.jar", "IntervalListTools"]
arguments: [valueFrom: "OUTPUT=$(runtime.outdir)/$(inputs.interval_list.nameroot).expanded.interval_list",
  "UNIQUE=TRUE"]
outputs:
  expanded_interval_list:
    type: File
    outputBinding:
      glob: "$(inputs.interval_list.nameroot).expanded.interval_list"
