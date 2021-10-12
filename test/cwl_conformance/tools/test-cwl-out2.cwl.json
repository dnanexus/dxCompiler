{
    "class": "CommandLineTool",
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
    "inputs": [],
    "outputs": [
        {
            "id": "#main/foo",
            "type": "File"
        }
    ],
    "arguments": [
        {
            "valueFrom": "echo foo > foo && echo '{\"foo\": {\"location\": \"file://$(runtime.outdir)/foo\", \"class\": \"File\"} }' > cwl.output.json\n",
            "shellQuote": false
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}