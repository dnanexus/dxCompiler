{
    "class": "CommandLineTool",
    "requirements": [
        {
            "listing": [
                "$(inputs.src)"
            ],
            "class": "InitialWorkDirRequirement"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "inputBinding": {
                "position": 1,
                "valueFrom": "$(self.nameroot).class"
            },
            "id": "#main/src"
        }
    ],
    "baseCommand": "touch",
    "outputs": [
        {
            "type": "File",
            "outputBinding": {
                "glob": "*.class"
            },
            "id": "#main/classfile"
        }
    ],
    "id": "#main",
    "cwlVersion": "v1.2"
}