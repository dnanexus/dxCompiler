{
    "class": "Workflow",
    "requirements": [
        {
            "class": "StepInputExpressionRequirement"
        },
        {
            "class": "MultipleInputFeatureRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [
        {
            "type": [
                "int",
                "string"
            ],
            "id": "#main/int_1"
        },
        {
            "type": [
                "int",
                "string"
            ],
            "id": "#main/int_2"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputSource": "#main/sum/result",
            "id": "#main/result"
        }
    ],
    "steps": [
        {
            "in": [
                {
                    "source": [
                        "#main/int_1",
                        "#main/int_2"
                    ],
                    "valueFrom": "${\n  var sum = 0;\n  for (var i = 0; i < self.length; i++){\n    sum += self[i];\n  };\n  return sum;\n}\n",
                    "id": "#main/sum/data"
                }
            ],
            "out": [
                "#main/sum/result"
            ],
            "run": {
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": "int",
                        "inputBinding": {},
                        "id": "#main/sum/run/data"
                    }
                ],
                "baseCommand": "echo",
                "stdout": "result.txt",
                "outputs": [
                    {
                        "type": "File",
                        "id": "#main/sum/run/result",
                        "outputBinding": {
                            "glob": "result.txt"
                        }
                    }
                ]
            },
            "id": "#main/sum"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}