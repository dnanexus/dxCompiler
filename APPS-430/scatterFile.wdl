version 1.0

task file_to_triple_string {
    input {
        File s
    }
    command <<<
        cat ~{s}
    >>>

    output {
        String res = read_string(stdout()) + read_string(stdout()) + read_string(stdout())
    }
}

task s_to_file {
    input {
        String s
    }

    command <<<
       echo ~{s} > ~{s}
    >>>

    output {
        File res = s
    }
}


workflow scatterFile {
    input {
        Array[String] strings
    }
    scatter (s in strings) {
        call s_to_file {
            input: s = s
        }
    }

    scatter (t in s_to_file.res) {
        call file_to_triple_string {
            input: s = t
        }
    }
    output {
        Array[String] tripled_string = file_to_triple_string.res
    }
}