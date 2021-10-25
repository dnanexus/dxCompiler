{
    "class": "Workflow",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "timelimit": 0,
            "class": "ToolTimeLimit"
        },
        {
            "enableReuse": false,
            "class": "WorkReuse"
        }
    ],
    "inputs": [
        {
            "type": [
                "null",
                "int"
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
                "class": "CommandLineTool",
                "baseCommand": [
                    "sleep",
                    "10"
                ],
                "inputs": [
                    {
                        "type": [
                            "null",
                            "int"
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
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}