version 1.0

task expressions_runtime_section {
    input {
    }
    command {
        echo "Something is happening"
    }
    runtime {
        docker: "dxcompilertest-there-is-no-such-image:123"
    }
    output {
    }
}
