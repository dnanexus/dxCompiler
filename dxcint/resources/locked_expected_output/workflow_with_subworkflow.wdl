version 1.1

import "subworkflow_with_defaults.wdl" as lib

workflow workflow_with_subworkflow {

    call increment_by_one { input: a = 10 }

    Int a = 3
    call lib.subworkflow_with_defaults as subworkflow_with_defaults {
        input: 
            a = a
            my_input=["hello", "world"]
    }

    output {
        Int r1 = increment_by_one.result
        Int r2 = subworkflow_with_defaults.result
        Array[String] array_output = subworkflow_with_defaults.array_output
    }
}

task increment_by_one {
    input {
        Int a
    }
    command {}
    output {
        Int result = a + 1
    }
}
