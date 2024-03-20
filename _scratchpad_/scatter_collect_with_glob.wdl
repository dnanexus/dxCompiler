version 1.1

workflow scatter_collect_with_glob {
  input {
    String input_string
    Int n = 8
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
    for i in {0..999}; do
      echo ~{input_string} > result_~{index}_${i}.txt
    done
  >>>
  output {
    Array[File] out = glob("result_*")
  }
}
