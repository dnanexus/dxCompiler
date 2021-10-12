{
    "class": "Workflow",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/my_file"
        }
    ],
    "steps": [
        {
            "run": {
                "class": "ExpressionTool",
                "requirements": [
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "inputs": [
                    {
                        "type": "File",
                        "loadContents": true,
                        "id": "#main/one/run/my_number"
                    }
                ],
                "outputs": [
                    {
                        "type": "int",
                        "id": "#main/one/run/my_int"
                    }
                ],
                "expression": "${ return { \"my_int\": parseInt(inputs.my_number.contents) }; }\n"
            },
            "in": [
                {
                    "source": "#main/my_file",
                    "id": "#main/one/my_number"
                }
            ],
            "out": [
                "#main/one/my_int"
            ],
            "id": "#main/one"
        }
    ],
    "outputs": [
        {
            "type": "int",
            "outputSource": "#main/one/my_int",
            "id": "#main/my_int"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}