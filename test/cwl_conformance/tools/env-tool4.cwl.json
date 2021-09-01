{
    "class": "CommandLineTool",
    "inputs": [
        {
            "type": "string",
            "id": "#main/in"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "out"
            },
            "id": "#main/out"
        }
    ],
    "requirements": [
        {
            "envDef": [
                {
                    "envValue": "conflict_original",
                    "envName": "TEST_ENV"
                }
            ],
            "class": "EnvVarRequirement"
        }
    ],
    "baseCommand": [
        "/bin/bash",
        "-c",
        "echo $TEST_ENV"
    ],
    "stdout": "out",
    "id": "#main",
    "cwlVersion": "v1.2"
}