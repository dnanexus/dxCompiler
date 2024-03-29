#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Picard MergeVcfs"
requirements:
- class: ResourceRequirement
  ramMin: 40000
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "broadinstitute/gatk:4.1.8.1"
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  merged_vcf_basename:
    type: string?
    default: "merged"
  sequence_dictionary:
    type:
    - string
    - File
    - 'null'
    inputBinding:
      position: 1
      prefix: "-D"
  vcfs:
    type:
      type: array
      items: File
      inputBinding:
        prefix: "-I"
    inputBinding:
      position: 2
baseCommand: ["/usr/bin/java", "-Xmx38g", "-jar", "/gatk/gatk.jar", "MergeVcfs"]
arguments: ["-O", "$(inputs.merged_vcf_basename).vcf.gz"]
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: $(inputs.merged_vcf_basename).vcf.gz
    secondaryFiles: [.tbi]
