{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "id": "#main/foo"
        }
    ],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "echo \"cow\" > \"$(runtime.outdir)/foo\" &&\necho \"moo\" > \"$(runtime.tmpdir)/foo\" &&\necho '{\"foo\": {\"path\": \"$(runtime.outdir)/foo\", \"class\": \"File\"} }' > cwl.output.json\n"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}