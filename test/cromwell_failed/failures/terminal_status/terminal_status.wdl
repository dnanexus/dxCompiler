task echo_sleep_exit {
    Int i
    command {
        echo ${i}
        sleep ${i}
        exit ${i}
    }
    runtime {
        continueOnReturnCode: [0, 2, 4]
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
        String o = read_string(stdout())
    }
}

# Shards 1 and 3 should fail, and 0, 2, 4 should succeed
# But all of them should reach a terminal status when the workflow ends
workflow terminal_status {
    scatter(i in range (5)) {
        call echo_sleep_exit { input: i = i }
    }
}
