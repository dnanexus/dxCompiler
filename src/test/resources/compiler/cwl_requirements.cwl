class: CommandLineTool
cwlVersion: v1.0

inputs:
  input: 'File[]?'

baseCommand: cat
arguments: [manifest.txt]

requirements:
  - class: DockerRequirement
    dockerPull: alpine
  - class: InitialWorkDirRequirement
    listing:
      - entryname: manifest.txt
        entry: |-
          ${
              var x = ""
              for (var i = 0; i < inputs.input.length; i++)  
              { 
                  x += inputs.input[i].path + "\n"; 
              }
              return x;
          }
        writable: false
  - class: InlineJavascriptRequirement

outputs:
  output: stdout

stdout: out.txt