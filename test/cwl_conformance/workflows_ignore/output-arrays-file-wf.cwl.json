{
    "class": "Workflow",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/i"
        }
    ],
    "outputs": [
        {
            "type": {
                "type": "array",
                "items": "File"
            },
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
                        "type": "File",
                        "loadContents": true,
                        "id": "#main/step1/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": {
                            "type": "array",
                            "items": "string"
                        },
                        "id": "#main/step1/run/o"
                    }
                ],
                "expression": "${return {'o': inputs.i.contents.split(\" \")};}\n"
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
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": {
                            "type": "array",
                            "items": "string"
                        },
                        "inputBinding": {
                            "position": 1
                        },
                        "id": "#main/step2/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": {
                            "type": "array",
                            "items": "File"
                        },
                        "outputBinding": {
                            "glob": "$(inputs.i)"
                        },
                        "id": "#main/step2/run/o"
                    }
                ],
                "baseCommand": "touch"
            },
            "id": "#main/step2"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}