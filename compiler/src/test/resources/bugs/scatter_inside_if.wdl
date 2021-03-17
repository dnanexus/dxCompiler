version 1.0

workflow test {
  input {
    Array[String]? str_input
  }

  if (defined(str_input)) {
    scatter (item in select_first([str_input])) {
      call show {
        input: a = item
      }
    }
    String abc = "def"
  }
}

task show {
  input {
    String? a
  }
  command {
    echo "~{a}"
  }
  output {}
}
