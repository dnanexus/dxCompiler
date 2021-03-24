{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "id": "#io-any-1.cwl/bar",
                    "type": "Any"
                }
            ],
            "outputs": [
                {
                    "id": "#io-any-1.cwl/t1",
                    "type": "Any",
                    "outputBinding": {
                        "outputEval": "$(inputs.bar.class || inputs.bar)"
                    }
                }
            ],
            "baseCommand": "true",
            "id": "#io-any-1.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "Any",
                    "id": "#main/bar"
                }
            ],
            "outputs": [
                {
                    "type": "Any",
                    "outputSource": "#main/step1/t1",
                    "id": "#main/t1"
                }
            ],
            "steps": [
                {
                    "in": [
                        {
                            "source": "#main/bar",
                            "id": "#main/step1/bar"
                        }
                    ],
                    "out": [
                        "#main/step1/t1"
                    ],
                    "run": "#io-any-1.cwl",
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.2"
}