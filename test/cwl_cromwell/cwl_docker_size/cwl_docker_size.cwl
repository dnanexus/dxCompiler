#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$graph:
- id: main
  class: CommandLineTool
  requirements:
  - class: InlineJavascriptRequirement
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
    NetworkAccess:
      networkAccess: true
    LoadListingRequirement:
      loadListing: deep_listing
  inputs:
  - id: INPUT
    type: File[]
  stdout: stdout
  outputs:
    stdout_output:
      type: string
      outputBinding:
        glob: stdout
        loadContents: true
        outputEval: $(self[0].contents.trim())
  arguments:
  - valueFrom: |
      ${
        if (inputs.INPUT.length == 0) {
          var cmd = ['echo', "no inputs"];
          return cmd
        }
        else {
          var cmd = ["echo", "execute"];
          var use_input = [];
          for (var i = 0; i < inputs.INPUT.length; i++) {
            var filesize = inputs.INPUT[i].size;
            use_input.push("order=".concat(filesize));
          }

          var run_cmd = cmd.concat(use_input);
          return run_cmd
        }

      }

  baseCommand: []
