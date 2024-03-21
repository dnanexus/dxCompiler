version 1.1

struct OutputStruct {
  File fileField
}

workflow scatter_collect_with_struct {
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
    echo ~{input_string} > struct_result_~{index}.txt
  >>>
  output {
    OutputStruct out = OutputStruct {
      fileField: "struct_result_~{index}.txt"
    }
  }
}
