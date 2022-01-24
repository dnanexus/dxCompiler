version 1.0

task division {
    input {
        Int number
    }
    command <<<
        echo ${number}
    >>>
    runtime {
        docker: "quay.io/ucsc_cgl/samtools"
    }
    output {
        Int result = 100/number
    }
}