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

    runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

workflow sizerelativepath {
    input {
        File file = "dx://file-GJ7VKG00yzZgfZzb1JgY2K2Q"
    }

    call print_size {
        input:
            file = file
    }
    output {
        String size_string = print_size.out
    }
}
