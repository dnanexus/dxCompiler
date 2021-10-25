{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "coresMin": "$(inputs.dir.listing[0].size)",
            "coresMax": "$(inputs.dir.listing[0].size)",
            "class": "ResourceRequirement"
        }
    ],
    "inputs": [
        {
            "type": "Directory",
            "id": "#main/dir"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "id": "#main/output",
            "outputBinding": {
                "glob": "cores.txt"
            }
        }
    ],
    "baseCommand": "echo",
    "stdout": "cores.txt",
    "arguments": [
        "$(runtime.cores)"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}