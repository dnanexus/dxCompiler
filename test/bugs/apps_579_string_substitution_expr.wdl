version 1.1

task apps_579_sub_t {
    input {
        String s = "aabb"
    }

    command <<<
    echo ~{sub(s, "bb$", "")}
    >>>

    output {
        String out = read_lines(stdout())[0]
    }
}

workflow apps_579_sub_wf {
    call apps_579_sub_t

    output {
        String out = apps_579_sub_t.out
    }
}
