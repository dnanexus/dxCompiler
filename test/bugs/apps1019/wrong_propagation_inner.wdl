version 1.0
task test_inner1 {
    input {
        String test_in
    }
    command <<<
        echo "~{test_in}" > "~{test_in}_inner_wf1.out2"
        exit 1
    >>>
    output {
        File test_out = "${test_in}_inner_wf1.out"
    }
}
task test_inner2 {
    input {
        File test_in
    }
    command <<<
        cat "~{test_in}" "~{test_in}" >  "~{test_in}_inner_wf2.out"
    >>>
    output {
        File test_out2 = "${test_in}_inner_wf2.out"
    }
}

workflow nested_inner {
    input {
        String s
    }
    call test_inner1 { input: test_in = s }
    call test_inner2 { input: test_in = test_inner1.test_out }
    output {
        File nested_inner_wf_out = test_inner2.test_out2
    }
}