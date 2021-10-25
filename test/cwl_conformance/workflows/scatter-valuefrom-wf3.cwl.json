{
    "$graph": [
        {
            "id": "#echo",
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#echo/echo_in1"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#echo/echo_in2"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#echo/first"
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
                    "id": "#echo/echo_out"
                }
            ],
            "baseCommand": "echo",
            "arguments": [
                "-n",
                "foo"
            ],
            "stdout": "step1_out"
        },
        {
            "id": "#main",
            "class": "Workflow",
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "name": "#main/inp1/instr",
                            "fields": [
                                {
                                    "name": "#main/inp1/instr/instr",
                                    "type": "string"
                                }
                            ]
                        }
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
                    "scatter": [
                        "#main/step1/echo_in1",
                        "#main/step1/echo_in2"
                    ],
                    "scatterMethod": "flat_crossproduct",
                    "in": [
                        {
                            "source": "#main/inp1",
                            "valueFrom": "$(self.instr)",
                            "id": "#main/step1/echo_in1"
                        },
                        {
                            "source": "#main/inp2",
                            "id": "#main/step1/echo_in2"
                        },
                        {
                            "source": "#main/inp1",
                            "valueFrom": "$(self[0].instr)",
                            "id": "#main/step1/first"
                        }
                    ],
                    "out": [
                        "#main/step1/echo_out"
                    ],
                    "run": "#echo",
                    "id": "#main/step1"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#main/step1/echo_out",
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "id": "#main/out"
                }
            ]
        }
    ],
    "cwlVersion": "v1.2"
}