{
    "class": "CommandLineTool",
    "doc": "Test of capturing stderr output in a docker container.",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "id": "#main/output_file",
            "outputBinding": {
                "glob": "std.err"
            }
        }
    ],
    "arguments": [
        {
            "valueFrom": "echo foo 1>&2",
            "shellQuote": false
        }
    ],
    "stderr": "std.err",
    "id": "#main",
    "cwlVersion": "v1.2"
}