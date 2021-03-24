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
            "type": "int",
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
                "class": "ExpressionTool",
                "inputs": [
                    {
                        "type": "int",
                        "id": "#main/sum/run/data"
                    }
                ],
                "outputs": [
                    {
                        "type": "int",
                        "id": "#main/sum/run/result"
                    }
                ],
                "expression": "${\n  return {\"result\": inputs.data};\n}\n"
            },
            "id": "#main/sum"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}