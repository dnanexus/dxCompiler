{
    "class": "Workflow",
    "doc": "Workflow without outputs.",
    "inputs": [
        {
            "type": "File",
            "id": "#main/file1"
        }
    ],
    "outputs": [],
    "steps": [
        {
            "in": [
                {
                    "source": "#main/file1",
                    "id": "#main/step0/file1"
                }
            ],
            "out": [],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "doc": "CommandLineTool without outputs.",
                "hints": [
                    {
                        "dockerPull": "debian:stretch-slim",
                        "class": "DockerRequirement"
                    }
                ],
                "inputs": [
                    {
                        "type": "File",
                        "label": "Input File",
                        "inputBinding": {
                            "position": 1
                        },
                        "id": "#main/step0/run/file1"
                    }
                ],
                "outputs": [],
                "baseCommand": "echo"
            },
            "id": "#main/step0"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}