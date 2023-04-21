version 1.0

task transubstantiation {

    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }

    command {
        echo 'foo' > output.txt
    }

    output {
        File out = "output.txt"
    }
}

workflow disproportionately {
    call transubstantiation

    output {
        File out = transubstantiation.out
    }
}

