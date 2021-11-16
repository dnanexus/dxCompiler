#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
inputs:
- id: pin
  type: string
  default: "v"
steps:
  threestep:
    run: cwl_three_step.cwl
    in:
    - id: pattern
      source: pin
    out: [wc-count]
outputs:
  count_output:
    type: int
    outputSource: threestep/wc-count
