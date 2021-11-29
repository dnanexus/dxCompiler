version 1.0

task s_double {
    input {
        String s
    }
    command <<<
        cd .
    >>>

    output {
        String res = s+s
    }
}


workflow reuseSomeFromScatter {
    input {
        Array[String] strings
    }
    scatter (s in strings) {
        call s_double {
            input: s = s
        }
    }
    output {
        Array[String] doubled_strings = s_double.res
    }
}