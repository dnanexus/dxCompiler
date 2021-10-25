version 1.0

workflow test {
  input {
    Array[File] file_input
  }
    scatter (item in select_first([file_input])) {
      call show {
        input: a = item
      }
    }
  output {
    Array[File]? result = show.output_file
  }
}

task show {
  input {
    File a
  }
  command {
    touch output_file.txt
    cat "~{a}" > output_file.txt
  }
  output {
    File output_file = "output_file.txt"
  }
}
