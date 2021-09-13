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
            "outputSource": "#main/step3/o",
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
                        "type": {
                            "type": "array",
                            "items": "int"
                        },
                        "id": "#main/step1/run/o"
                    }
                ],
                "expression": "${return {'o': Array.apply(null, {length: inputs.i}).map(Number.call, Number)};}\n"
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
                        "type": {
                            "type": "array",
                            "items": "int"
                        },
                        "id": "#main/step2/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": {
                            "type": "array",
                            "items": "int"
                        },
                        "id": "#main/step2/run/o"
                    }
                ],
                "expression": "${return {'o': inputs.i.map(function(x) { return (x + 1) * 2; })};}\n"
            },
            "id": "#main/step2"
        },
        {
            "in": [
                {
                    "source": "#main/step2/o",
                    "id": "#main/step3/i"
                }
            ],
            "out": [
                "#main/step3/o"
            ],
            "run": {
                "class": "ExpressionTool",
                "inputs": [
                    {
                        "type": {
                            "type": "array",
                            "items": "int"
                        },
                        "id": "#main/step3/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": "int",
                        "id": "#main/step3/run/o"
                    }
                ],
                "expression": "${return {'o': inputs.i.reduce(function(a, b) { return a + b; })};}\n"
            },
            "id": "#main/step3"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}