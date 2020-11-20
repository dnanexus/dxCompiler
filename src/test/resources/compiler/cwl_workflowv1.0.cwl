class: Workflow
cwlVersion: v1.0
id: my_workflow

inputs:
  one: string

steps:
  first:
    run: first.cwl
    in:
      { id: input1, source: one }
    out: [result]
  second:
    run: second.cwl
    in: {id: input2: source: first/result }
    out: []
outputs: []