{
    "class": "CommandLineTool",
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "out-filelist.txt",
                    "entry": "${\n  var ls = \"\";\n  for (var i = 1; i < 10000; i++) {\n    ls += \"example_input_file\"+i+\".txt\\n\";\n  }\n  return ls;\n}"
                }
            ],
            "class": "InitialWorkDirRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [],
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "out-filelist.txt"
            },
            "id": "#main/filelist"
        }
    ],
    "baseCommand": "true",
    "id": "#main",
    "cwlVersion": "v1.2"
}