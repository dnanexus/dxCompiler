{
    "class": "CommandLineTool",
    "requirements": [
        {
            "listing": [
                {
                    "entry": "$(inputs.file1)",
                    "entryname": "bob.txt"
                }
            ],
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "#main/file1"
        }
    ],
    "outputs": [],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "test \"$(inputs.file1.path)\" = \"$(runtime.outdir)/bob.txt\"\n"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}