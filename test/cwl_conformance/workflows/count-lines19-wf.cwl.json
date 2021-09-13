{
    "$graph": [
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "File",
                    "id": "#main/file1"
                }
            ],
            "outputs": [
                {
                    "type": "int",
                    "outputSource": "#main/step1/output",
                    "id": "#main/count_output"
                }
            ],
            "steps": [
                {
                    "run": "#wc3-tool.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/file1"
                            ],
                            "linkMerge": "merge_nested",
                            "id": "#main/step1/file1"
                        }
                    ],
                    "out": [
                        "#main/step1/output"
                    ],
                    "id": "#main/step1"
                }
            ],
            "id": "#main"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {},
                    "id": "#wc3-tool.cwl/file1"
                }
            ],
            "outputs": [
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "output.txt",
                        "loadContents": true,
                        "outputEval": "${\n  var s = self[0].contents.split(/\\r?\\n/);\n  return parseInt(s[s.length-2]);\n}\n"
                    },
                    "id": "#wc3-tool.cwl/output"
                }
            ],
            "stdout": "output.txt",
            "baseCommand": "wc",
            "id": "#wc3-tool.cwl"
        }
    ],
    "cwlVersion": "v1.2"
}