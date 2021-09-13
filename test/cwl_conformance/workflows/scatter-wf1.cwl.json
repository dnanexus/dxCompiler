{
    "class": "Workflow",
    "inputs": [
        {
            "type": {
                "type": "array",
                "items": "string"
            },
            "id": "#main/inp"
        }
    ],
    "outputs": [
        {
            "type": {
                "type": "array",
                "items": "string"
            },
            "outputSource": "#main/step1/echo_out",
            "id": "#main/out"
        }
    ],
    "requirements": [
        {
            "class": "ScatterFeatureRequirement"
        }
    ],
    "steps": [
        {
            "in": [
                {
                    "source": "#main/inp",
                    "id": "#main/step1/echo_in"
                }
            ],
            "out": [
                "#main/step1/echo_out"
            ],
            "scatter": "#main/step1/echo_in",
            "run": {
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": "string",
                        "inputBinding": {},
                        "id": "#main/step1/run/echo_in"
                    }
                ],
                "outputs": [
                    {
                        "type": "string",
                        "outputBinding": {
                            "glob": "step1_out",
                            "loadContents": true,
                            "outputEval": "$(self[0].contents)"
                        },
                        "id": "#main/step1/run/echo_out"
                    }
                ],
                "baseCommand": "echo",
                "arguments": [
                    "-n",
                    "foo"
                ],
                "stdout": "step1_out"
            },
            "id": "#main/step1"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}