version 1.1

task apps_579_boolean_t {
    input {
        Boolean bool_flag = true
    }

    command <<<
    if_true="x"
    if_false="y"
    echo ~{true="$if_true" false="$if_false" bool_flag}
    >>>

    output {
        String out = read_lines(stdout())[0]
    }
}

workflow apps_579_boolean_wf {
    call apps_579_boolean_t

    output {
        String out = apps_579_boolean_t.out
    }
}
