{
    "$graph": [
        {
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "string",
                    "id": "#env-tool2.cwl/in"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "out"
                    },
                    "id": "#env-tool2.cwl/out"
                }
            ],
            "hints": [
                {
                    "envDef": [
                        {
                            "envValue": "$(inputs.in)",
                            "envName": "TEST_ENV"
                        }
                    ],
                    "class": "EnvVarRequirement"
                }
            ],
            "baseCommand": [
                "/bin/sh",
                "-c",
                "echo $TEST_ENV"
            ],
            "stdout": "out",
            "id": "#env-tool2.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "string",
                    "id": "#main/in"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/step1/out",
                    "id": "#main/out"
                }
            ],
            "requirements": [
                {
                    "envDef": [
                        {
                            "envValue": "override",
                            "envName": "TEST_ENV"
                        }
                    ],
                    "class": "EnvVarRequirement"
                }
            ],
            "steps": [
                {
                    "run": "#env-tool2.cwl",
                    "in": [
                        {
                            "source": "#main/in",
                            "id": "#main/step1/in"
                        }
                    ],
                    "out": [
                        "#main/step1/out"
                    ],
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.2"
}