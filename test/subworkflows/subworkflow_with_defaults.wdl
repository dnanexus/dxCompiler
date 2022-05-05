version 1.1

# A workflow that has a required argument with a default.
# There was a bug (https://github.com/dnanexus/dxWDL/issues/284) where
# a call failed because it did not supplied [b].
workflow subworkflow_with_defaults {
    input {
        Int a
        Int b = 10
        Array[String] my_input = []
    }

    call print_array {
        input: a = my_input
    }

    output {
        Int result = a + b
        Array[String] array_output = my_input
        String print_result = print_array.result
    }
}

# Prints the values of a String array
task print_array {
     input {
         Array[String] a
     }

     command {
         echo ${sep=' ' a}
     }

     output {
         String result = read_string(stdout())
     }
}
