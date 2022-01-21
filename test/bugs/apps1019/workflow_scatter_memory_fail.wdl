version 1.0

import "memoryfail.wdl" as mem
workflow workflow_scatter_memory_fail {
    input {
        Array[Int] numbers
    }

    scatter (i in numbers) {
        call mem.memoryfail {
            input: number = i
        }
    }
    output {
        Array[File] result = memoryfail.result
    }
}