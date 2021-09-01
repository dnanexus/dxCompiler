{
    "class": "CommandLineTool",
    "id": "#main",
    "baseCommand": [
        "ls"
    ],
    "inputs": [
        {
            "id": "#input_list",
            "type": {
                "type": "array",
                "items": "Directory"
            }
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
                "glob": [
                    "testdir/a",
                    "rec/B"
                ]
            }
        }
    ],
    "label": "stage-array-dirs.cwl",
    "requirements": [
        {
            "class": "InitialWorkDirRequirement",
            "listing": [
                "$(inputs.input_list)"
            ]
        }
    ],
    "cwlVersion": "v1.2"
}