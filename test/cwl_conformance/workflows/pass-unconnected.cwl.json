{
    "$graph": [
        {
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "Any",
                    "inputBinding": {},
                    "id": "#echo-tool.cwl/in"
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
                    "id": "#echo-tool.cwl/out"
                }
            ],
            "baseCommand": "echo",
            "stdout": "out.txt",
            "id": "#echo-tool.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "string",
                    "default": "hello inp1",
                    "id": "#main/inp1"
                },
                {
                    "type": "string",
                    "default": "hello inp2",
                    "id": "#main/inp2"
                }
            ],
            "outputs": [
                {
                    "type": "string",
                    "outputSource": "#main/step1/out",
                    "id": "#main/out"
                }
            ],
            "steps": [
                {
                    "in": [
                        {
                            "source": "#main/inp1",
                            "id": "#main/step1/in"
                        },
                        {
                            "source": "#main/inp2",
                            "id": "#main/step1/in2"
                        }
                    ],
                    "out": [
                        "#main/step1/out"
                    ],
                    "run": "#echo-tool.cwl",
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.2"
}