version 1.0

# Tests: scatter, write_tsv, inline python code, ragged string arrays
#

# Generate sets of intervals for scatter-gathering over chromosomes
task createTsv {
    # Use python to create a string parsed into a wdl Array[Array[String]]
    command<<<
    python3 <<CODE
    tsv_list = []
    ll = [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
    for l in ll:
      tsv_list.append('\t'.join(l))
    tsv_string = '\n'.join(tsv_list)
    print(tsv_string)
    CODE
    >>>

    output {
      Array[Array[String]] result = read_tsv(stdout())
    }
}

task processTsv {
    input {
        Array[Array[String]] words
    }
    Array[Array[String]] first_line = [words[0]]

    command {
        cat ${write_tsv(first_line)}
    }
    output {
      # It used to be true that read_string only read the first
      # line of the file, but that is not actually correct according
      # to the spec. read_string now returns the entire file as a single
      # line, with only the trailing newline stripped off.
      String result = read_string(stdout())
    }
}

task processLine {
    input {
        Array[String] line
    }
    command {
        echo ${sep=' INPUT=' line}
    }
    output {
        String result = read_string(stdout())
    }
}

task writeLines {
    input {
        String x
        String y
    }

    command <<<
    cat ~{write_lines(["~{x} ~{y}"])}
    >>>

    output {
        String result = read_lines(stdout())[0]
    }

    runtime {
        docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
}

workflow strings {
    input {
    }

    # Ragged array of strings
    call createTsv
    call processTsv { input: words = createTsv.result }
    scatter (x in createTsv.result) {
        call processLine {input : line=x}
    }
    call writeLines {
        input: x = "hello", y = "buddy"
    }
    output {
        String result1 = processTsv.result
        Array[String] result2 = processLine.result
        String result3 = writeLines.result
    }
}
