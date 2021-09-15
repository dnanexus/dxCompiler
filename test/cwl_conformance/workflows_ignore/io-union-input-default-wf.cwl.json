{
    "class": "Workflow",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [
        {
            "type": [
                "File",
                "null",
                "string"
            ],
            "default": "the default value",
            "id": "#main/bar"
        }
    ],
    "outputs": [
        {
            "type": "string",
            "outputSource": "#main/step1/o",
            "id": "#main/o"
        }
    ],
    "steps": [
        {
            "in": [
                {
                    "source": "#main/bar",
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
                        "type": [
                            "File",
                            "null",
                            "string"
                        ],
                        "id": "#main/step1/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": "string",
                        "id": "#main/step1/run/o"
                    }
                ],
                "expression": "${return {'o': (inputs.i.class || inputs.i)};}\n"
            },
            "id": "#main/step1"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}