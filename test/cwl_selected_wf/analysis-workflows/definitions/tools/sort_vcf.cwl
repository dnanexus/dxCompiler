#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Sort VCF"
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
  vcf:
    type: File
    inputBinding:
      prefix: "I="
  reference_dict:
    type: File?
    inputBinding:
      prefix: "SEQUENCE_DICTIONARY="
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "SortVcf"]
arguments: ["O=", valueFrom: $(runtime.outdir)/sorted.vcf]
outputs:
  sorted_vcf:
    type: File
    outputBinding:
      glob: "sorted.vcf"
