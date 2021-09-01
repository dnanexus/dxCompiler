{
    "class": "CommandLineTool",
    "doc": "Test of capturing stderr output.",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "error.txt"
            },
            "id": "#main/output_file"
        }
    ],
    "arguments": [
        {
            "valueFrom": "echo foo 1>&2",
            "shellQuote": false
        }
    ],
    "stderr": "error.txt",
    "id": "#main",
    "cwlVersion": "v1.2"
}