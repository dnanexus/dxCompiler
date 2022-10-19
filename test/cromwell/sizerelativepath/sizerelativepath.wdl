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
        File file = "dx://file-GJ74pxQ0yzZy59jJ4p5FY930"
    }

    call print_size {
        input:
            file = file
    }
    output {
        String size_string = print_size.out
    }
}
