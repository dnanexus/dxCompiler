#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Custom allele frequency filter"
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "mgibio/vep_helper-cwl:vep_101.0_v1"
- class: ResourceRequirement
  ramMin: 4000
- class: StepInputExpressionRequirement
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
  maximum_population_allele_frequency:
    type: float
    inputBinding:
      valueFrom: |
        ${
            return [
                "--filter",
                [
                    inputs.field_name, "<", inputs.maximum_population_allele_frequency,
                    "or not", inputs.field_name
                ].join(" ")
            ]
        }
      position: 2
  field_name:
    type: string
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
arguments: [valueFrom: $(inputs.vcf.path), valueFrom: $(runtime.outdir)/annotated.af_filtered.vcf,
  "/usr/bin/perl", "/opt/vep/src/ensembl-vep/filter_vep", "--format", "vcf", "-o",
  valueFrom: $(runtime.outdir)/annotated.af_filtered.vcf]
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: "annotated.af_filtered.vcf"
