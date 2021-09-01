{
    "class": "CommandLineTool",
    "doc": "Must fail if entryname starts with ../\n",
    "inputs": [
        {
            "type": "File",
            "id": "#main/filelist"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "head.txt"
            },
            "id": "#main/head"
        }
    ],
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "../input/stuff.txt",
                    "entry": "$(inputs.filelist)"
                }
            ],
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "head -n10 ../input/stuff.txt > head.txt"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}