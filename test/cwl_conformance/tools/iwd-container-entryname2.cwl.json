{
    "class": "CommandLineTool",
    "doc": "Must fail if entryname is an absolute path and DockerRequirement is\nnot in the 'requirements' section.\n",
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
                    "entryname": "/tmp2j3y7rpb/input/stuff.txt",
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
            "valueFrom": "head -n10 /tmp2j3y7rpb/input/stuff.txt > /output/head.txt"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}