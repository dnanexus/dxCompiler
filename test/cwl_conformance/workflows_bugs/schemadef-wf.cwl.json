{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "SchemaDefRequirement",
                    "types": [
                        {
                            "name": "#schemadef-type.yml/HelloType",
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#schemadef-type.yml/HelloType/a",
                                    "type": "string"
                                },
                                {
                                    "name": "#schemadef-type.yml/HelloType/b",
                                    "type": "string"
                                }
                            ]
                        }
                    ],
                    "id": "#schemadef-type.yml",
                    "name": "#schemadef-type.yml"
                }
            ],
            "inputs": [
                {
                    "id": "#schemadef-tool.cwl/hello",
                    "type": "#schemadef-type.yml/HelloType",
                    "inputBinding": {
                        "valueFrom": "$(self.a)/$(self.b)"
                    }
                }
            ],
            "outputs": [
                {
                    "id": "#schemadef-tool.cwl/output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "output.txt"
                    }
                }
            ],
            "stdout": "output.txt",
            "baseCommand": "echo",
            "id": "#schemadef-tool.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "$import": "#schemadef-type.yml"
                }
            ],
            "inputs": [
                {
                    "type": "#schemadef-type.yml/HelloType",
                    "id": "#main/hello"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/step1/output",
                    "id": "#main/output"
                }
            ],
            "steps": [
                {
                    "in": [
                        {
                            "source": "#main/hello",
                            "id": "#main/step1/hello"
                        }
                    ],
                    "out": [
                        "#main/step1/output"
                    ],
                    "run": "#schemadef-tool.cwl",
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.2"
}