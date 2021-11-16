task invalid_return_code {
    command <<<
        echo successful
    >>>
    output {
        String successful = read_string(stdout())
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
        continueOnReturnCode: 1
    }
}

workflow invalid_return_code_wf {
    call invalid_return_code
}
