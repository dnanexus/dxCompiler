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
            "id": "#main/i"
        }
    ],
    "outputs": [
        {
            "type": "int",
            "outputSource": "#main/step1/o",
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
                "expression": "${return {'o': inputs.i * 2};}\n"
            },
            "id": "#main/step1"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}