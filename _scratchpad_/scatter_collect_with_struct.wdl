version 1.1

struct OutputStruct {
  File fileField
}

workflow scatter_collect_no_glob {
  input {
    String input_string
    Int n = 3
  }
  scatter (i in range(n)) {
    call scatter_task {
      input:
        input_string = input_string,
        index = i
    }
  }
}

task scatter_task {
  input {
    String input_string
    Int index
  }
  command <<<
    echo ~{input_string} > no_glob_result_~{index}.txt
  >>>
  output {
    OutputStruct out = OutputStruct {
      fileField: "no_glob_result_~{index}.txt"
    }
  }
}
