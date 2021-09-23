{
    "requirements": [
        {
            "class": "DockerRequirement",
            "dockerPull": "debian:stretch-slim"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "class": "CommandLineTool",
    "inputs": [
        {
            "type": {
                "type": "array",
                "items": "File"
            },
            "secondaryFiles": [
                {
                    "pattern": ".fai",
                    "required": true
                },
                {
                    "pattern": ".crai",
                    "required": false
                },
                {
                    "pattern": ".bai",
                    "required": false
                },
                {
                    "pattern": "${ if (inputs.require_dat) {return '.dat'} else {return null} }",
                    "required": null
                },
                {
                    "pattern": "${ return null; }",
                    "required": null
                },
                {
                    "pattern": ".dat2",
                    "required": "$(inputs.require_dat)"
                }
            ],
            "id": "#main/fasta_path"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "id": "#main/require_dat"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "fai.list"
            },
            "secondaryFiles": [
                {
                    "pattern": ".bai",
                    "required": false
                },
                {
                    "pattern": "${ return null }"
                }
            ],
            "id": "#main/bai_list"
        }
    ],
    "arguments": [
        {
            "valueFrom": "${ var fai_list = \"\"; for (var i = 0; i < inputs.fasta_path.length; i ++) { fai_list += \" cat \" + inputs.fasta_path[i].path +\".fai\" + \" >> fai.list && \" } return fai_list.slice(0,-3) }",
            "position": 1,
            "shellQuote": false
        }
    ],
    "baseCommand": [],
    "id": "#main",
    "cwlVersion": "v1.2"
}