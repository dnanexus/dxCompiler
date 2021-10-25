{
    "class": "CommandLineTool",
    "requirements": [
        {
            "dockerPull": "debian:stretch-slim",
            "dockerOutputDirectory": "/other",
            "class": "DockerRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "thing"
            },
            "id": "#main/thing"
        }
    ],
    "baseCommand": [
        "touch",
        "/other/thing"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}