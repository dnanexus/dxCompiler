{
    "class": "Workflow",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [
        {
            "type": "int",
            "default": 4,
            "id": "#main/i"
        }
    ],
    "outputs": [
        {
            "type": "int",
            "outputSource": "#main/step2/o",
            "id": "#main/o"
        }
    ],
    "steps": [
        {
            "in": [
                {
                    "source": "#main/i",
                    "id": "#main/step1/i"
                }
            ],
            "out": [
                "#main/step1/o"
            ],
            "run": {
                "class": "ExpressionTool",
                "inputs": [
                    {
                        "type": "int",
                        "id": "#main/step1/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": "int",
                        "id": "#main/step1/run/o"
                    }
                ],
                "expression": "${return {'o': (inputs.i || 2)};}\n"
            },
            "id": "#main/step1"
        },
        {
            "in": [
                {
                    "source": "#main/step1/o",
                    "id": "#main/step2/i"
                }
            ],
            "out": [
                "#main/step2/o"
            ],
            "run": {
                "class": "ExpressionTool",
                "inputs": [
                    {
                        "type": "int",
                        "id": "#main/step2/run/i"
                    },
                    {
                        "type": "int",
                        "default": 5,
                        "id": "#main/step2/run/i2"
                    }
                ],
                "outputs": [
                    {
                        "type": "int",
                        "id": "#main/step2/run/o"
                    }
                ],
                "expression": "${return {'o': inputs.i * 2 + inputs.i2};}\n"
            },
            "id": "#main/step2"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}