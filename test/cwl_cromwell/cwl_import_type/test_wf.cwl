#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
requirements:
- class: SchemaDefRequirement
  types:
  - $import: capture_kit.yml
inputs:
- id: bam
  type: string
- id: capture_kit
  type: capture_kit.yml#capture_kit
steps:
- id: touch_bam
  run: touch.cwl
  in:
  - id: input
    source: bam
  out:
  - id: output
outputs:
- id: output_bam
  type: File
  outputSource: touch_bam/output
