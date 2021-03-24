{
    "class": "Workflow",
    "doc": "Workflow without inputs.",
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "outputSource": "#main/step0/output",
            "id": "#main/output"
        }
    ],
    "steps": [
        {
            "in": [],
            "out": [
                "#main/step0/output"
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "doc": "CommandLineTool without inputs.",
                "hints": [
                    {
                        "dockerPull": "debian:stretch-slim",
                        "class": "DockerRequirement"
                    }
                ],
                "inputs": [],
                "outputs": [
                    {
                        "type": "File",
                        "outputBinding": {
                            "glob": "output"
                        },
                        "id": "#main/step0/run/output"
                    }
                ],
                "baseCommand": [
                    "echo",
                    "cwl"
                ],
                "stdout": "output"
            },
            "id": "#main/step0"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}