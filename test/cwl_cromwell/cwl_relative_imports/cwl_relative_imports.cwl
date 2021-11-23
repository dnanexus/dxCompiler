#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  class: Workflow
  inputs: []
  outputs:
    letters:
      type: string
      outputSource: globSort/letters
  steps:
    globSort:
      in: []
      run: cwl_relative_imports_glob_sort.cwl
      out: [letters]
