#!/usr/bin/env cwl-runner

cwlVersion: v1.2

requirements:
  - class: DockerRequirement
    dockerPull: ubuntu:focal-20200423
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement

class: CommandLineTool

inputs:
  - id: create
    type: boolean
    default: true
    inputBinding:
      prefix: --create
      position: 0

  - id: dirname
    type: string

  - id: gzip
    type: boolean
    default: false
    inputBinding:
      prefix: --gzip
      position: 2

  - id: transform
    type: ["null", string]
      
  - id: input
    type:
      type: array
      items: File

arguments:
  - valueFrom: $(inputs.dirname).tar
    prefix: --file
    position: 1

  - valueFrom: |
      ${
        var inp_str = "";
        for (var i = 0; i < inputs.input.length; i++) {
          var cmd = " -C " + inputs.input[i].dirname + " " + inputs.input[i].basename;
          inp_str = inp_str.concat(cmd);
         }
        return inp_str; 
       }
    position: 99
    shellQuote: false

  - valueFrom: |
      ${
          var cmd = "--transform 's,^," + inputs.dirname + "/,'";
          return cmd;
      }
    position: 90
    shellQuote: false
    
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.dirname).tar

baseCommand: [tar]
