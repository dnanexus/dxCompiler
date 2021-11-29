version 1.0

task identity2 {
    input {
        File file
    }
    command <<<
        cd .
    >>>

    output {
        File res = file
    }
}

task slice_file2 {
    input {
        Int bytes
    }

    command <<<
        head -c ~{bytes} </dev/urandom > random.out
        split -b1 -a6 --additional-suffix=.bam random.out outf
    >>>

    output {
        Array[File] slices = glob("outf*")
    }
}


workflow slice_and_identity2 {
    input {
        Int number_of_files
        File some_file
    }

    call slice_file2 {
        input: bytes = number_of_files
    }

    scatter (slice in slice_file2.slices) {
        call identity2 {
            input: file = slice,
        }
    }

    output {
        Array[File] counts = identity2.res
    }
}