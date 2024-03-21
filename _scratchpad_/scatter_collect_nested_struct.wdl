version 1.1

struct InnerStruct {
  File fileField
}

struct OuterStruct {
  File fileField
  InnerStruct structField
}

workflow scatter_collect_nested_struct {
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
    echo ~{input_string} > nested_array_result_~{index}_1.txt
    echo ~{input_string} > nested_array_result_~{index}_2.txt
  >>>
  output {
    OuterStruct out = OuterStruct {
      fileField: "struct_result_~{index}_1.txt",
      structField: InnerStruct {
        fileField: "struct_result_~{index}_2.txt"
      }
    }
  }
}
