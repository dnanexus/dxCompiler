#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Coding Variant filter"
requirements:
- class: DockerRequirement
  dockerPull: "mgibio/vep_helper-cwl:vep_101.0_v1"
- class: ResourceRequirement
  ramMin: 4000
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  vcf:
    type: File
    inputBinding:
      prefix: "-i"
      position: 1
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
arguments: [valueFrom: $(inputs.vcf.path), valueFrom: $(runtime.outdir)/annotated.coding_variant_filtered.vcf,
  "/usr/bin/perl", "/opt/vep/src/ensembl-vep/filter_vep", "--format", "vcf", "-o",
  valueFrom: $(runtime.outdir)/annotated.coding_variant_filtered.vcf, "--ontology",
  "--filter", "Consequence is coding_sequence_variant"]
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: "annotated.coding_variant_filtered.vcf"
