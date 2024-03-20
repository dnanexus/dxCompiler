version 1.1

workflow scatter_collect_with_glob {
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
    for i in {0..99}; do
      echo ~{input_string} > glob_result_~{index}_${i}.txt
    done
  >>>
  output {
    Array[File] out = glob("glob_result_*")
  }
}
