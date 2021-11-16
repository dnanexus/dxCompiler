task invalid_runtime_attributes {
    command { # NOOP }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
        continueOnReturnCode: "oops"
    }
}

workflow invalid_runtime_attributes_wf {
    call invalid_runtime_attributes
}
