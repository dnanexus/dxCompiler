version 1.0

import "workflow_scatter_division.wdl" as wf_div
workflow workflow_calling_workflow {
    input {
        Array[Array[Int]] s
    }
    scatter (i in s) {
        call wf_div.workflow_scatter_division {
            input: numbers = i
        }
    }
    output {
        Array[Array[Int]] out = workflow_scatter_division.result
    }
}