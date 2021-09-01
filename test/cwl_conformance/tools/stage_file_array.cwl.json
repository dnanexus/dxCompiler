{
    "class": "CommandLineTool",
    "id": "#main",
    "label": "Stage File Array",
    "arguments": [
        "ls"
    ],
    "inputs": [
        {
            "id": "#input_list",
            "type": {
                "type": "array",
                "items": "File"
            },
            "secondaryFiles": [
                {
                    "pattern": ".sec",
                    "required": null
                }
            ]
        }
    ],
    "outputs": [
        {
            "id": "#output",
            "type": {
                "type": "array",
                "items": "File"
            },
            "outputBinding": {
                "glob": "input_dir/*"
            }
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "InitialWorkDirRequirement",
            "listing": [
                {
                    "entryname": "input_dir",
                    "entry": "${ return {class: 'Directory', listing: inputs.input_list} }"
                }
            ]
        }
    ],
    "cwlVersion": "v1.2"
}