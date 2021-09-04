{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/inf"
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
        "$(inputs.inf.dirname)/xtestdir",
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