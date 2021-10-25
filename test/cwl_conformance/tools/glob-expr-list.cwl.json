{
    "class": "CommandLineTool",
    "inputs": [
        {
            "type": {
                "type": "array",
                "items": "string"
            },
            "inputBinding": {
                "position": 1
            },
            "id": "#main/ids"
        }
    ],
    "outputs": [
        {
            "type": {
                "type": "array",
                "items": "File"
            },
            "outputBinding": {
                "glob": "$(inputs.ids)"
            },
            "id": "#main/files"
        }
    ],
    "baseCommand": "touch",
    "id": "#main",
    "cwlVersion": "v1.2"
}