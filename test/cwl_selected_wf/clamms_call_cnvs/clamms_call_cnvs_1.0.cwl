baseCommand: []
class: CommandLineTool
cwlVersion: v1.0
id: clamms_call_cnvs
inputs:
  cnv_rate:
    default: 1.0e-07
    doc: Reflects the expected rate of CNVs across the exome. Lower values increase
      sensitivity.
    inputBinding:
      position: 4
      prefix: --cnv_rate
    type: double?
  reference_model_file:
    doc: Pre-computed model file fit on a reference panel
    inputBinding:
      position: 2
      prefix: --reference_model_file
    type: File
  sample_normalized_coverage_bed_file:
    doc: BED file of pre-computed normalized coverage for the test sample
    inputBinding:
      position: 1
      prefix: --sample_normalized_coverage_bed_file
    type: File
  sex:
    doc: '"M" or "F" for assigning expected ploidy on sex chromosomes. If not provided,
      CNVs will not be called on sex chormosomes.'
    inputBinding:
      position: 3
      prefix: --sex
    type: string?
label: 'CLAMMS: Call CNVs for a sample using a pre-computed reference model'
outputs:
  cnv_bed_file:
    doc: BED file of CNVs called by CLAMMS on the test sample
    outputBinding:
      glob: cnv_bed_file/*
    type: File
requirements:
- class: DockerRequirement
  dockerOutputDirectory: /data/out
  dockerPull: pfda2dockstore/clamms_call_cnvs:12
s:author:
  class: s:Person
  s:name: Evan Maxwell