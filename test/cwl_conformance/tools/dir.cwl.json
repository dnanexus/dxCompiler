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
        "cd",
        "$(inputs.indir.path)",
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