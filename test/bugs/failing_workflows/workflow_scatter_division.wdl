version 1.0

import "division.wdl" as div
workflow workflow_scatter_division {
    input {
        Array[Int] numbers
    }
    scatter (i in numbers) {
        call div.division {
            input: number = i
        }
    }
    output {
        Array[Int] result = division.result
    }
}