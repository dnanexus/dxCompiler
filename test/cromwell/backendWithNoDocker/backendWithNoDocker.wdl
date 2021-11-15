task dockerhub {
    command {
        echo "bonjour tout le monde !"
    }
    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
        backend: "LocalNoDocker"
    }
}

workflow backend_with_no_docker {
    call dockerhub
}
