task diff2files {
    File a
    File b

    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    command {
        diff ${a} ${b} | wc -l
    }
    output {
        Int result = read_int(stdout())
    }
    parameter_meta {
        a : "stream"
        b : "stream"
    }
}
