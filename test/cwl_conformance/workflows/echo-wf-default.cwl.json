{
    "$graph": [
        {
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "string",
                    "default": "tool_default",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#echo-tool-default.cwl/in"
                }
            ],
            "outputs": [
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "out.txt",
                        "loadContents": true,
                        "outputEval": "$(self[0].contents)"
                    },
                    "id": "#echo-tool-default.cwl/out"
                }
            ],
            "baseCommand": [
                "echo",
                "-n"
            ],
            "stdout": "out.txt",
            "id": "#echo-tool-default.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [],
            "steps": [
                {
                    "run": "#echo-tool-default.cwl",
                    "in": [
                        {
                            "default": "workflow_default",
                            "id": "#main/step1/in"
                        }
                    ],
                    "out": [
                        "#main/step1/out"
                    ],
                    "id": "#main/step1"
                }
            ],
            "outputs": [
                {
                    "type": "string",
                    "outputSource": "#main/step1/out",
                    "id": "#main/default_output"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.2"
}