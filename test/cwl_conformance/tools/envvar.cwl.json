{
    "class": "CommandLineTool",
    "inputs": [],
    "outputs": [],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "arguments": [
        "echo",
        {
            "valueFrom": "\"HOME=$HOME\"",
            "shellQuote": false
        },
        {
            "valueFrom": "\"TMPDIR=$TMPDIR\"",
            "shellQuote": false
        },
        {
            "valueFrom": "&&",
            "shellQuote": false
        },
        "test",
        {
            "valueFrom": "\"$HOME\"",
            "shellQuote": false
        },
        "=",
        "$(runtime.outdir)",
        "-a",
        {
            "valueFrom": "\"$TMPDIR\"",
            "shellQuote": false
        },
        "=",
        "$(runtime.tmpdir)"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}