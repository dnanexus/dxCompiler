{
    "class": "CommandLineTool",
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "results"
            },
            "id": "#main/results"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "hints": [
        {
            "dockerPull": "debian:stretch-slim",
            "class": "DockerRequirement"
        }
    ],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "echo HOME=$HOME TMPDIR=$TMPDIR > log\nif [ \"$HOME\" = \"$(runtime.outdir)\" ] && [ \"$TMPDIR\" = \"$(runtime.tmpdir)\" ]\nthen\n    echo success > results\nelse\n    echo failure > results\nfi\n"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}