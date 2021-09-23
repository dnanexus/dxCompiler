{
    "class": "CommandLineTool",
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "emptyWritableDir",
                    "writable": true,
                    "entry": "$({class: 'Directory', listing: []})"
                }
            ],
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "Directory",
            "outputBinding": {
                "glob": "emptyWritableDir"
            },
            "id": "#main/out"
        }
    ],
    "arguments": [
        "touch",
        "emptyWritableDir/blurg"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}