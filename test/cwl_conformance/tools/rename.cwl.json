{
    "class": "CommandLineTool",
    "baseCommand": "true",
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "$(inputs.newname)",
                    "entry": "$(inputs.srcfile)"
                }
            ],
            "class": "InitialWorkDirRequirement"
        }
    ],
    "inputs": [
        {
            "type": "string",
            "id": "#main/newname"
        },
        {
            "type": "File",
            "id": "#main/srcfile"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.newname)"
            },
            "id": "#main/outfile"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}