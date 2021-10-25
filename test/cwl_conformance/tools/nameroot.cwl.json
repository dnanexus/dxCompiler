{
    "class": "CommandLineTool",
    "inputs": [
        {
            "type": "File",
            "id": "#main/file1"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "id": "#main/b",
            "outputBinding": {
                "glob": "$(inputs.file1.nameroot).xtx"
            }
        }
    ],
    "stdout": "$(inputs.file1.nameroot).xtx",
    "baseCommand": [],
    "arguments": [
        "echo",
        "$(inputs.file1.basename)",
        "$(inputs.file1.nameroot)",
        "$(inputs.file1.nameext)"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}