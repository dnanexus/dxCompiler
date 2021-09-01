{
    "class": "CommandLineTool",
    "inputs": [
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
                        "name": "#main/record_input/f1"
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
                        "name": "#main/record_input/f2"
                    }
                ]
            },
            "id": "#main/record_input"
        }
    ],
    "outputs": [],
    "baseCommand": "test",
    "arguments": [
        "-f",
        "$(inputs.record_input.f1.path).s2",
        "-a",
        "-f",
        "$(inputs.record_input.f2[0].path).s3",
        "-a",
        "-f",
        "$(inputs.record_input.f2[1].path).s3"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}