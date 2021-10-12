{
    "requirements": [
        {
            "class": "DockerRequirement",
            "dockerPull": "debian:stretch-slim"
        },
        {
            "class": "InitialWorkDirRequirement",
            "listing": [
                "$(inputs.INPUT)"
            ]
        }
    ],
    "class": "CommandLineTool",
    "inputs": [
        {
            "id": "#main/INPUT",
            "type": "File"
        }
    ],
    "outputs": [
        {
            "id": "#main/OUTPUT",
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.INPUT.basename)"
            },
            "secondaryFiles": [
                {
                    "pattern": ".fai",
                    "required": null
                }
            ]
        }
    ],
    "arguments": [
        {
            "valueFrom": "$(inputs.INPUT.basename).fai",
            "position": 0
        }
    ],
    "baseCommand": [
        "touch"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}