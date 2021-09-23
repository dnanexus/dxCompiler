{
    "class": "CommandLineTool",
    "doc": "Create a file under adir/, symlink it to working directory (./) and glob symlink. The executor should resolve this symlink",
    "hints": [
        {
            "dockerPull": "alpine",
            "class": "DockerRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "symlink.txt"
            },
            "id": "#main/output_file"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "arguments": [
        "mkdir",
        "adir",
        {
            "valueFrom": " && ",
            "shellQuote": false
        },
        "echo",
        "Who's gonna drive you home",
        {
            "valueFrom": "> adir/original.txt",
            "shellQuote": false
        },
        {
            "valueFrom": " && ",
            "shellQuote": false
        },
        "ln",
        "-s",
        "adir/original.txt",
        "symlink.txt"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}