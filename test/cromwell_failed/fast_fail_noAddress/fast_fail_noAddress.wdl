task task_with_noAddress {

    command { echo "hello, world" }

    output { String s = read_string(stdout()) }

    runtime {
      docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
      # Our network is not configured specially, so this should cause this task to fail almost immediately.
      noAddress: true
    }
}

workflow fast_fail_noAddress {
    call task_with_noAddress
}
