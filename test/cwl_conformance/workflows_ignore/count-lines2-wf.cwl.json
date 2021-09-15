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
            "id": "#main/file1"
        }
    ],
    "outputs": [
        {
            "type": "int",
            "outputSource": "#main/step2/parseInt_output",
            "id": "#main/count_output"
        }
    ],
    "steps": [
        {
            "in": [
                {
                    "source": "#main/file1",
                    "id": "#main/step1/wc_file1"
                }
            ],
            "out": [
                "#main/step1/wc_output"
            ],
            "run": {
                "id": "#main/step1/run/wc",
                "class": "CommandLineTool",
                "inputs": [
                    {
                        "type": "File",
                        "inputBinding": {},
                        "id": "#main/step1/run/wc/wc_file1"
                    }
                ],
                "outputs": [
                    {
                        "type": "File",
                        "outputBinding": {
                            "glob": "output.txt"
                        },
                        "id": "#main/step1/run/wc/wc_output"
                    }
                ],
                "stdout": "output.txt",
                "baseCommand": "wc"
            },
            "id": "#main/step1"
        },
        {
            "in": [
                {
                    "source": "#main/step1/wc_output",
                    "id": "#main/step2/parseInt_file1"
                }
            ],
            "out": [
                "#main/step2/parseInt_output"
            ],
            "run": {
                "class": "ExpressionTool",
                "inputs": [
                    {
                        "type": "File",
                        "loadContents": true,
                        "id": "#main/step2/run/parseInt_file1"
                    }
                ],
                "outputs": [
                    {
                        "type": "int",
                        "id": "#main/step2/run/parseInt_output"
                    }
                ],
                "expression": "${return {'parseInt_output': parseInt(inputs.parseInt_file1.contents)};}\n"
            },
            "id": "#main/step2"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}