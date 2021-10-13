{
    "$graph": [
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "File",
                    "id": "#count-lines1-wf-noET.cwl/file1"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#count-lines1-wf-noET.cwl/step1/output",
                    "id": "#count-lines1-wf-noET.cwl/wc_output"
                }
            ],
            "steps": [
                {
                    "run": "#wc-tool.cwl",
                    "in": [
                        {
                            "source": "#count-lines1-wf-noET.cwl/file1",
                            "id": "#count-lines1-wf-noET.cwl/step1/file1"
                        }
                    ],
                    "out": [
                        "#count-lines1-wf-noET.cwl/step1/output"
                    ],
                    "id": "#count-lines1-wf-noET.cwl/step1"
                }
            ],
            "id": "#count-lines1-wf-noET.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "File",
                    "id": "#main/file1"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/step1/wc_output",
                    "id": "#main/wc_output"
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "steps": [
                {
                    "run": "#count-lines1-wf-noET.cwl",
                    "in": [
                        {
                            "source": "#main/file1",
                            "id": "#main/step1/file1"
                        }
                    ],
                    "out": [
                        "#main/step1/wc_output"
                    ],
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        },
        {
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "File",
                    "id": "#wc-tool.cwl/file1"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "output"
                    },
                    "id": "#wc-tool.cwl/output"
                }
            ],
            "baseCommand": [
                "wc",
                "-l"
            ],
            "stdin": "$(inputs.file1.path)",
            "stdout": "output",
            "id": "#wc-tool.cwl"
        }
    ],
    "cwlVersion": "v1.2"
}