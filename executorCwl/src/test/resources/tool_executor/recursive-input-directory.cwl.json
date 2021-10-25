{
    "class": "CommandLineTool",
    "requirements": [
        {
            "listing": [
                {
                    "entry": "$(inputs.input_dir)",
                    "entryname": "work_dir",
                    "writable": true
                }
            ],
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "ShellCommandRequirement"
        }
    ],
    "stdout": "output.txt",
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "touch work_dir/e;\nif [ ! -w work_dir ]; then echo work_dir not writable; fi;\nif [ -L work_dir ]; then echo work_dir is a symlink; fi;\nif [ ! -w work_dir/a ]; then echo work_dir/a not writable; fi;\nif [ -L work_dir/a ]; then echo work_dir/a is a symlink; fi;\nif [ ! -w work_dir/c ]; then echo work_dir/c not writable; fi;\nif [ -L work_dir/c ]; then echo work_dir/c is a symlink; fi;\nif [ ! -w work_dir/c/d ]; then echo work_dir/c/d not writable; fi;\nif [ -L work_dir/c/d ]; then echo work_dir/c/d is a symlink; fi;\nif [ ! -w work_dir/e ]; then echo work_dir/e not writable; fi;\nif [ -L work_dir/e ]; then echo work_dir/e is a symlink ; fi;\n"
        }
    ],
    "inputs": [
        {
            "type": "Directory",
            "id": "#main/input_dir"
        }
    ],
    "outputs": [
        {
            "type": "Directory",
            "outputBinding": {
                "glob": "work_dir"
            },
            "id": "#main/output_dir"
        },
        {
            "type": "File",
            "id": "#main/test_result",
            "outputBinding": {
                "glob": "output.txt"
            }
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}