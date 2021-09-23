{
    "class": "CommandLineTool",
    "inputs": [],
    "outputs": [
        {
            "type": {
                "type": "record",
                "fields": [
                    {
                        "type": "File",
                        "secondaryFiles": {
                            "pattern": ".s2",
                            "required": null
                        },
                        "outputBinding": {
                            "glob": "A"
                        },
                        "name": "#main/record_output/f1"
                    },
                    {
                        "type": {
                            "type": "array",
                            "items": "File"
                        },
                        "secondaryFiles": {
                            "pattern": ".s3",
                            "required": null
                        },
                        "outputBinding": {
                            "glob": [
                                "B",
                                "C"
                            ]
                        },
                        "name": "#main/record_output/f2"
                    }
                ]
            },
            "id": "#main/record_output"
        }
    ],
    "baseCommand": "touch",
    "arguments": [
        "A",
        "A.s2",
        "B",
        "B.s3",
        "C",
        "C.s3"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}