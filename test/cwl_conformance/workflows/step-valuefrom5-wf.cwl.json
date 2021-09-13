{
    "class": "Workflow",
    "requirements": [
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/file1"
        }
    ],
    "outputs": [
        {
            "type": "string",
            "outputSource": "#main/step1/echo_out",
            "id": "#main/val1"
        },
        {
            "type": "string",
            "outputSource": "#main/step2/echo_out",
            "id": "#main/val2"
        }
    ],
    "steps": [
        {
            "run": {
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": "string",
                        "inputBinding": {},
                        "id": "#main/step1/run/name"
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
                    },
                    {
                        "type": "File",
                        "outputBinding": {
                            "glob": "step1_out"
                        },
                        "id": "#main/step1/run/echo_out_file"
                    }
                ],
                "baseCommand": "echo",
                "stdout": "step1_out"
            },
            "in": [
                {
                    "source": "#main/file1",
                    "valueFrom": "$(self.basename)",
                    "id": "#main/step1/name"
                }
            ],
            "out": [
                "#main/step1/echo_out",
                "#main/step1/echo_out_file"
            ],
            "id": "#main/step1"
        },
        {
            "run": {
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": "string",
                        "inputBinding": {},
                        "id": "#main/step2/run/name"
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
                        "id": "#main/step2/run/echo_out"
                    }
                ],
                "baseCommand": "echo",
                "stdout": "step1_out"
            },
            "in": [
                {
                    "source": "#main/step1/echo_out_file",
                    "valueFrom": "$(self.basename)",
                    "id": "#main/step2/name"
                }
            ],
            "out": [
                "#main/step2/echo_out"
            ],
            "id": "#main/step2"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}