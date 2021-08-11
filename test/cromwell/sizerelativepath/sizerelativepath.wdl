version 1.0

task print_size {
    input {
        File file
    }
    Int bytes = ceil(size(file))

    command {
        echo ~{bytes}
    }

    output {
        String out = read_string(stdout())
    }

    runtime {docker: "ubuntu:latest"}
}

workflow sizerelativepath {
    input {
        File file = "dx://file-G3f1Y5j0yzZy5G13GVY8J78Q"
    }

    call print_size {
        input:
            file = file
    }
    output {
        String size_string = print_size.out
    }
}
