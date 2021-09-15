{
    "class": "Workflow",
    "requirements": [
        {
            "coresMin": 4,
            "coresMax": 4,
            "class": "ResourceRequirement"
        }
    ],
    "inputs": [],
    "steps": [
        {
            "requirements": [
                {
                    "coresMin": 1,
                    "coresMax": 1,
                    "class": "ResourceRequirement"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "inputs": [],
                "outputs": [
                    {
                        "type": "File",
                        "id": "#main/step1/run/output",
                        "outputBinding": {
                            "glob": "cores.txt"
                        }
                    }
                ],
                "baseCommand": "echo",
                "stdout": "cores.txt",
                "arguments": [
                    "$(runtime.cores)"
                ]
            },
            "in": [],
            "out": [
                "#main/step1/output"
            ],
            "id": "#main/step1"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputSource": "#main/step1/output",
            "id": "#main/out"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}