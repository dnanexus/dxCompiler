{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [
        {
            "type": {
                "name": "#main/irec/irec",
                "type": "record",
                "fields": [
                    {
                        "name": "#main/irec/irec/ifoo",
                        "type": "File",
                        "inputBinding": {
                            "position": 2
                        }
                    },
                    {
                        "name": "#main/irec/irec/ibar",
                        "type": "File",
                        "inputBinding": {
                            "position": 6
                        }
                    }
                ]
            },
            "id": "#main/irec"
        }
    ],
    "outputs": [
        {
            "type": {
                "name": "#main/orec/orec",
                "type": "record",
                "fields": [
                    {
                        "name": "#main/orec/orec/ofoo",
                        "type": "File",
                        "outputBinding": {
                            "glob": "foo"
                        }
                    },
                    {
                        "name": "#main/orec/orec/obar",
                        "type": "File",
                        "outputBinding": {
                            "glob": "bar"
                        }
                    }
                ]
            },
            "id": "#main/orec"
        }
    ],
    "arguments": [
        {
            "valueFrom": "cat",
            "position": 1
        },
        {
            "valueFrom": "> foo",
            "position": 3,
            "shellQuote": false
        },
        {
            "valueFrom": "&&",
            "position": 4,
            "shellQuote": false
        },
        {
            "valueFrom": "cat",
            "position": 5
        },
        {
            "valueFrom": "> bar",
            "position": 7,
            "shellQuote": false
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}