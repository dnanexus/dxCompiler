{
    "class": "CommandLineTool",
    "requirements": [
        {
            "dockerPull": "debian:stable-slim",
            "class": "DockerRequirement"
        },
        {
            "listing": "${\n   return [{\"class\": \"Directory\",\n            \"basename\": \"subdir\",\n            \"listing\": [ inputs.filelist ]\n            }]}\n",
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
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
                "glob": "subdir/$(inputs.filelist.basename)"
            },
            "id": "#main/same"
        }
    ],
    "baseCommand": "echo",
    "id": "#main",
    "cwlVersion": "v1.2"
}