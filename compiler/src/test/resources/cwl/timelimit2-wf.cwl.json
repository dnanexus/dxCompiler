{
    "class": "Workflow",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "timelimit": 5,
            "class": "ToolTimeLimit"
        }
    ],
    "inputs": [
        {
            "type": [
                "null",
                "string"
            ],
            "id": "#main/i"
        }
    ],
    "outputs": [
        {
            "type": [
                "null",
                "string"
            ],
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
                "class": "CommandLineTool",
                "baseCommand": [
                    "sleep",
                    "3"
                ],
                "inputs": [
                    {
                        "type": [
                            "null",
                            "string"
                        ],
                        "id": "#main/step1/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": [
                            "null",
                            "string"
                        ],
                        "outputBinding": {
                            "outputEval": "$(\"time passed\")"
                        },
                        "id": "#main/step1/run/o"
                    }
                ]
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
                "baseCommand": [
                    "sleep",
                    "3"
                ],
                "inputs": [
                    {
                        "type": [
                            "null",
                            "string"
                        ],
                        "id": "#main/step2/run/i"
                    }
                ],
                "outputs": [
                    {
                        "type": [
                            "null",
                            "string"
                        ],
                        "outputBinding": {
                            "outputEval": "$(\"time passed\")"
                        },
                        "id": "#main/step2/run/o"
                    }
                ]
            },
            "id": "#main/step2"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}