#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "collect alignment summary metrics"
requirements:
- class: ResourceRequirement
  ramMin: 18000
- class: DockerRequirement
  dockerPull: "broadinstitute/picard:2.23.6"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  bam:
    type: File
    inputBinding:
      prefix: "INPUT="
    secondaryFiles: [^.bai]
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "REFERENCE_SEQUENCE="
  metric_accumulation_level:
    type: string
    inputBinding:
      prefix: "METRIC_ACCUMULATION_LEVEL="
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectAlignmentSummaryMetrics"]
arguments: ["OUTPUT=", valueFrom: $(runtime.outdir)/$(inputs.bam.nameroot).AlignmentSummaryMetrics.txt]
outputs:
  alignment_summary_metrics:
    type: File
    outputBinding:
      glob: "$(inputs.bam.nameroot).AlignmentSummaryMetrics.txt"
