version 1.0

task dyno_runtime_reuse {
    input {
        Int i
    }
    command <<<
        echo "~{i}" > ~{i}
        cat /etc/os-release >>~{i}
    >>>
    output {
        File dyno_out="${i}"
    }
    runtime {
        docker: if i > 2 then "ubuntu:focal" else "ubuntu:bionic"
        memory: if i > 2 then "32GB" else "8GB"
    }
}