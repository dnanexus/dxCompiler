version 1.1

workflow scatter_collect_nested_array {
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
    echo ~{input_string} > nested_array_result_~{index}_3.txt
    echo ~{input_string} > nested_array_result_~{index}_4.txt
  >>>
  output {
    Array[Array[File]] out = [
      ["nested_array_result_~{index}_1.txt", "nested_array_result_~{index}_2.txt"],
      ["nested_array_result_~{index}_3.txt", "nested_array_result_~{index}_4.txt"]
    ]
  }
}
