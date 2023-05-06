task apps_912_t01 {
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


workflow apps_912 {
    input {
        Array[String] strings
    }
    scatter (s in strings) {
        call apps_912_t01 {
            input: s = s
        }
    }
    output {
        Array[String] doubled_strings = apps_912_t01.res
    }
}