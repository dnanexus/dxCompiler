{
    "class": "Workflow",
    "requirements": [
        {
            "class": "StepInputExpressionRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "MultipleInputFeatureRequirement"
        }
    ],
    "inputs": [
        {
            "type": "int",
            "id": "#main/a"
        },
        {
            "type": "int",
            "id": "#main/b"
        }
    ],
    "outputs": [
        {
            "type": "string",
            "outputSource": "#main/step1/echo_out",
            "id": "#main/val"
        }
    ],
    "steps": [
        {
            "run": {
                "id": "#main/step1/run/echo",
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": "int",
                        "inputBinding": {},
                        "id": "#main/step1/run/echo/c"
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
                        "id": "#main/step1/run/echo/echo_out"
                    }
                ],
                "baseCommand": "echo",
                "stdout": "step1_out"
            },
            "in": [
                {
                    "source": [
                        "#main/a",
                        "#main/b"
                    ],
                    "valueFrom": "$(self[0] + self[1])",
                    "id": "#main/step1/c"
                }
            ],
            "out": [
                "#main/step1/echo_out"
            ],
            "id": "#main/step1"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}