task empty_glob {
    command {
        echo "hello"
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
    output {
        Array[File] empty_glob = glob("*.txt")
        Int nbGlobbed = length(empty_glob)
    }
}

workflow wdl_empty_glob {
    call empty_glob
    output {
        Int wf_out = empty_glob.nbGlobbed
    }
}
