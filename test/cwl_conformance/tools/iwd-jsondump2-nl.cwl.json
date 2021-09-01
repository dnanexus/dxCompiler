{
    "class": "CommandLineTool",
    "requirements": [
        {
            "listing": [
                {
                    "entryname": "out-filelist.json",
                    "entry": "${\n  var ls = [];\n  for (var i = 1; i < 10000; i++) {\n    ls.push(\"example_input_file\"+i+\".txt\");\n  }\n  return {\"filelist\": ls};\n}\n"
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
                "glob": "out-filelist.json"
            },
            "id": "#main/filelist"
        }
    ],
    "arguments": [
        "true"
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}