#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: relative_imports
  class: Workflow
  inputs: []
  outputs:
    letters:
      type: string
      outputSource: globSort/letters
  steps:
    globSort:
      in: []
      run: cwl_glob_sort.cwl
      out: [letters]