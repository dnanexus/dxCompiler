version 1.0

task identity {
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

task slice_file {
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


workflow slice_and_identity {
    input {
        Int number_of_files
    }

    call slice_file {
        input: bytes = number_of_files
    }

    scatter (slice in slice_file.slices) {
        call identity {
            input: file = slice,
        }
    }

    output {
        Array[File] counts = identity.res
    }
}