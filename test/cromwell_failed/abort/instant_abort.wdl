task aborted {
    command {
        echo "Instant abort"
    }
    
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
}

workflow instant_abort {
    call aborted
}
