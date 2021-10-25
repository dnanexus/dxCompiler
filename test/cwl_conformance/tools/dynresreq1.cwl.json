{
    "class": "CommandLineTool",
    "requirements": [
        {
            "coresMin": "$(inputs.special_file.size)",
            "coresMax": "$(inputs.special_file.size)",
            "class": "ResourceRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/special_file"
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