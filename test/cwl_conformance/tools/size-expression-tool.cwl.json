{
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "class": "CommandLineTool",
    "inputs": [
        {
            "id": "#main/input",
            "type": {
                "type": "array",
                "items": "File"
            }
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "output.txt"
            },
            "id": "#main/output_file"
        }
    ],
    "arguments": [
        {
            "valueFrom": "${\n  var cmd = [\"echo\"];\n  if (inputs.input.length == 0) {\n     cmd.push('no_inputs');\n  }\n  else {\n    for (var i = 0; i < inputs.input.length; i++) {\n      var filesize = inputs.input[i].size;\n      if (filesize == 0) {\n        cmd.push(\"empty_file\");\n      } else if (filesize <= 16) {\n        cmd.push(\"small_file\");\n      } else {\n        cmd.push(\"big_file\")\n      }\n    }\n  }\n  return cmd;\n}\n"
        }
    ],
    "baseCommand": [],
    "stdout": "output.txt",
    "id": "#main",
    "cwlVersion": "v1.2"
}