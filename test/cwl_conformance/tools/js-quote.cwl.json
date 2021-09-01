{
    "class": "CommandLineTool",
    "label": "Check for a JS quoting bug",
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "file.txt",
                    "entry": "${return 'quote \"' + inputs.quote + '\"'}\n"
                },
                {
                    "entryname": "script.sh",
                    "entry": "set -xe\ncat file.txt\n"
                }
            ],
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "coresMin": 2,
            "ramMin": 1000,
            "class": "ResourceRequirement"
        }
    ],
    "inputs": [
        {
            "type": "string",
            "default": "Hello",
            "id": "#main/quote"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "file.txt"
            },
            "id": "#main/out"
        }
    ],
    "baseCommand": [
        "echo"
    ],
    "arguments": [],
    "id": "#main",
    "cwlVersion": "v1.2"
}