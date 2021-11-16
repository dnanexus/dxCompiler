# A task which takes an input, allows this to be called using a WDL function in the input block
task wdl_function_locations_i {
    Array[Int] array
    command {
      echo ${sep=":" array} > out
    }
    runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }
    output { String out = read_string("out") }
}

task wdl_function_locations_ii {
    Int x
    command { 
      echo ${x} > out
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
        # A (pure) WDL function in the runtime attributes:
        continueOnReturnCode: range(3)
    }
    output { String out = read_string("out") }
}

workflow wdl_function_locations {
    call wdl_function_locations_i { input: array = range(3) }
    scatter (x in range(2)) {
        call wdl_function_locations_ii { input: x = x }
    }
}
