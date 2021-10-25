{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [
        {
            "type": "Directory",
            "inputBinding": {
                "prefix": "cd",
                "position": -1
            },
            "id": "#main/indir"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "output.txt"
            },
            "id": "#main/outlist"
        }
    ],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "&&"
        },
        "find",
        ".",
        {
            "shellQuote": false,
            "valueFrom": "|"
        },
        "sort"
    ],
    "stdout": "output.txt",
    "id": "#main",
    "cwlVersion": "v1.2"
}