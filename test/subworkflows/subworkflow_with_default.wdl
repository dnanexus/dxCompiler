version 1.0

import "Foo.wdl" as lib

workflow subworkflow_with_default {

    call Inc { input: a = 10 }

    Int a = 3
    call lib.Foo as Foo {
        input: 
            a = a
            my_input=["hello"]
    }

    output {
        Int r1 = Inc.result
        Int r2 = Foo.result
        Array[String] my_output = Foo.my_output
    }
}

task Inc {
    input {
        Int a
    }
    command {}
    output {
        Int result = a + 1
    }
}
