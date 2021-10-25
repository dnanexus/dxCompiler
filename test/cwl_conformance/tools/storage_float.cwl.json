{
    "class": "CommandLineTool",
    "requirements": [
        {
            "ramMin": 254.1,
            "ramMax": 254.9,
            "tmpdirMin": 255.1,
            "tmpdirMax": 255.9,
            "outdirMin": 256.1,
            "outdirMax": 256.9,
            "class": "ResourceRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "id": "#main/output",
            "outputBinding": {
                "glob": "values.txt"
            }
        }
    ],
    "baseCommand": "echo",
    "stdout": "values.txt",
    "arguments": [
        "$(runtime.ram)",
        "$(runtime.tmpdirSize)",
        "$(runtime.outdirSize)"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}