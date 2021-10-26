#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: EnvVarRequirement
  envDef:
    TMPDIR: /tmp
hints:
  DockerRequirement:
    dockerPull: karchinlab/opencravat
  NetworkAcess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  input:
    type: File
    inputBinding:
      position: 1
      shellQuote: false
  modulesDirectory:
    type: Directory
    inputBinding:
      prefix: '--system-option modules_dir='
      separate: false
      position: 2
      shellQuote: false
  genome:
    type:
      type: enum
      symbols:
      - hg38
      - hg19
    inputBinding:
      prefix: -l
      position: 3
      shellQuote: false
    default: hg38
baseCommand: ['oc', 'run']
arguments:
- prefix: -d
  valueFrom: '.'
  shellQuote: false
- prefix: --endat
  valueFrom: postaggregator
  shellQuote: false
outputs:
  db:
    type: File
    outputBinding:
      glob: $(inputs.input.basename).sqlite
  log:
    type: File
    outputBinding:
      glob: $(inputs.input.basename).log
  err:
    type: File
    outputBinding:
      glob: $(inputs.input.basename).err
