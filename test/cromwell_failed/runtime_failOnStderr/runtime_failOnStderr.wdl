task failOnStderr {
    command <<<
        echo "OH NO!" >&2
    >>>
    output {
        String ohno = read_string(stderr())
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
        failOnStderr: true
    }
}

workflow runtime_failOnStderr {
    call failOnStderr
}
