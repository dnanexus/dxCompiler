{
    "class": "Workflow",
    "inputs": [
        {
            "type": {
                "type": "array",
                "items": {
                    "type": "record",
                    "name": "#main/inp/instr",
                    "fields": [
                        {
                            "name": "#main/inp/instr/instr",
                            "type": "string"
                        }
                    ]
                }
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
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "steps": [
        {
            "in": [
                {
                    "source": "#main/inp",
                    "valueFrom": "$(self.instr)",
                    "id": "#main/step1/echo_in"
                },
                {
                    "source": "#main/inp",
                    "valueFrom": "$(inputs.echo_in.instr)",
                    "id": "#main/step1/first"
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
                        "inputBinding": {
                            "position": 2
                        },
                        "id": "#main/step1/run/echo_in"
                    },
                    {
                        "type": "string",
                        "inputBinding": {
                            "position": 1
                        },
                        "id": "#main/step1/run/first"
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