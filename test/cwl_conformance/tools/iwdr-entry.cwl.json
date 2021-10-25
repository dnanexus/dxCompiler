{
    "class": "CommandLineTool",
    "baseCommand": [
        "cat",
        "example.conf"
    ],
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "example.conf",
                    "entry": "CONFIGVAR=$(inputs.message)\n"
                }
            ],
            "class": "InitialWorkDirRequirement"
        }
    ],
    "inputs": [
        {
            "type": "string",
            "id": "#main/message"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "example.conf"
            },
            "id": "#main/out"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}