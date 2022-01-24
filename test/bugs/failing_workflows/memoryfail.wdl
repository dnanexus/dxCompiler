version 1.0

# Task will cause out of memory error if called with argument == 0
task memoryfail {
    input {
        Int number
    }
    command <<<
        echo "running task with.... ~{number}"
        if [[ "~{number}" -eq 0 ]]; then
            fork() {
                fork | fork &
            }
        fork
        fi
        touch file.out
        echo $number > file.out
    >>>
    runtime {
        docker: "quay.io/ucsc_cgl/samtools"
    }
    output {
        File result = "file.out"
    }
}
