#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "pindel somatic filter v1"
requirements:
- class: ResourceRequirement
  ramMin: 16000
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: mgibio/cle:v1.3.1
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
  pindel_output_summary:
    type: File
arguments: ["/usr/bin/perl", "/usr/bin/write_pindel_filter_config.pl", $(inputs.pindel_output_summary.path),
  $(inputs.reference), $(runtime.outdir), {valueFrom: " && ", shellQuote: false},
  "/usr/bin/perl", "/usr/bin/somatic_indelfilter.pl", "filter.config"]
outputs:
  vcf:
    type: File
    outputBinding:
      glob: "pindel.out.vcf"
