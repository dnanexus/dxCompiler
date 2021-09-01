{
    "class": "CommandLineTool",
    "doc": "When executing in a container, entryname can have an absolute path\nto a mount location inside the container.\n",
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
            "dockerPull": "debian:10",
            "dockerOutputDirectory": "/output",
            "class": "DockerRequirement"
        },
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