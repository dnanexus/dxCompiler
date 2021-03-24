{
    "$graph": [
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/file1"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "outputSource": "#main/step1/output",
                    "id": "#main/count_output"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "steps": [
                {
                    "run": "#wc2-tool.cwl",
                    "scatter": "#main/step1/file1",
                    "in": [
                        {
                            "source": "#main/file1",
                            "id": "#main/step1/file1"
                        }
                    ],
                    "out": [
                        "#main/step1/output"
                    ],
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "id": "#wc2-tool.cwl/file1",
                    "type": "File",
                    "inputBinding": {}
                }
            ],
            "outputs": [
                {
                    "id": "#wc2-tool.cwl/output",
                    "type": "int",
                    "outputBinding": {
                        "glob": "output.txt",
                        "loadContents": true,
                        "outputEval": "$(parseInt(self[0].contents))"
                    }
                }
            ],
            "stdout": "output.txt",
            "baseCommand": "wc",
            "id": "#wc2-tool.cwl"
        }
    ],
    "cwlVersion": "v1.2"
}