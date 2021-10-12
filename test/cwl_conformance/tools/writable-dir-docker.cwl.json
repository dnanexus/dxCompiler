{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "InitialWorkDirRequirement",
            "listing": [
                {
                    "entryname": "emptyWritableDir",
                    "entry": "$({class: 'Directory', listing: []})",
                    "writable": true
                }
            ]
        }
    ],
    "hints": [
        {
            "class": "DockerRequirement",
            "dockerPull": "alpine"
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