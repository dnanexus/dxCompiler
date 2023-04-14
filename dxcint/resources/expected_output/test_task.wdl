version 1.0


task test_task {
    input {
        Int integer2 = 3332
        Int integer1 = 3333
        String string = "test_string3"
        Boolean boolean = true
    }

    command {
        echo ~{integer2}
        echo ~{integer1}
        echo ~{string}
        echo ~{boolean}
    }


    output {
        Int integer1_res = integer1
        Int integer2_res = integer2
        String string_res = string
        Boolean boolean_res = boolean
    }
}
