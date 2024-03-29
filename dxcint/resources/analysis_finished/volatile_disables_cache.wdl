version 1.0

workflow volatile_disables_cache {
    call volatile_task
    output {
        Int random = volatile_task.out
    }
}

task volatile_task {
    meta {
        description: "Do something stochastic. We don't want this to call cache!"
        volatile: true
    }
    command <<<
        echo $RANDOM
    >>>
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
        Int out = read_int(stdout())
    }
}
