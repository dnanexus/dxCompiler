{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "InitialWorkDirRequirement",
            "listing": "$(inputs.indir.listing)"
        }
    ],
    "inputs": [
        {
            "type": "Directory",
            "loadListing": "shallow_listing",
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
        "find",
        "-L",
        ".",
        "!",
        "-path",
        "*.txt",
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