{
    "class": "CommandLineTool",
    "hints": [
        {
            "class": "DockerRequirement",
            "dockerPull": "python:2-slim"
        }
    ],
    "requirements": [
        {
            "listing": [
                {
                    "entry": "$(inputs.infile)",
                    "entryname": "bob.txt",
                    "writable": true
                }
            ],
            "class": "InitialWorkDirRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/infile"
        }
    ],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "bob.txt"
            },
            "id": "#main/outfile"
        }
    ],
    "baseCommand": "python",
    "arguments": [
        "-c",
        "f = open(\"bob.txt\", \"r+\")\nf.seek(8)\nf.write(\"Bob.    \")\nf.close()\n"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}