{
    "class": "Workflow",
    "requirements": [
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "inputs": [],
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
                        "type": "string",
                        "inputBinding": {},
                        "id": "#main/step1/run/echo/a"
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
                    "valueFrom": "moocow",
                    "id": "#main/step1/a"
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