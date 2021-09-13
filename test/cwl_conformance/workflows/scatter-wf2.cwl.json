{
    "class": "Workflow",
    "inputs": [
        {
            "type": {
                "type": "array",
                "items": "string"
            },
            "id": "#main/inp1"
        },
        {
            "type": {
                "type": "array",
                "items": "string"
            },
            "id": "#main/inp2"
        }
    ],
    "outputs": [
        {
            "outputSource": "#main/step1/echo_out",
            "type": {
                "type": "array",
                "items": {
                    "type": "array",
                    "items": "string"
                }
            },
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
                    "source": "#main/inp1",
                    "id": "#main/step1/echo_in1"
                },
                {
                    "source": "#main/inp2",
                    "id": "#main/step1/echo_in2"
                }
            ],
            "out": [
                "#main/step1/echo_out"
            ],
            "scatter": [
                "#main/step1/echo_in1",
                "#main/step1/echo_in2"
            ],
            "scatterMethod": "nested_crossproduct",
            "run": {
                "class": "CommandLineTool",
                "id": "#main/step1/run/step1command",
                "inputs": [
                    {
                        "type": "string",
                        "inputBinding": {},
                        "id": "#main/step1/run/step1command/echo_in1"
                    },
                    {
                        "type": "string",
                        "inputBinding": {},
                        "id": "#main/step1/run/step1command/echo_in2"
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
                        "id": "#main/step1/run/step1command/echo_out"
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