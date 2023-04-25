task file_outputs_from_input_t01 {
  String outputName1
  String outputName2
  String outputName3

  command {
    echo "foo" > ${outputName1}
    echo "bar" > ${outputName2}
    echo "baz" > ${outputName3}.txt
  }
  output {
    File foo = outputName1
    File bar = "${outputName2}"
    File baz = "${outputName3}.txt"
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow file_outputs_from_input {
  call file_outputs_from_input_t01 { input: outputName1 = "___something___", outputName2 = "___anything___", outputName3 = "___nothing___" }
  output {
    String foo = read_string(file_outputs_from_input_t01.foo)
    String bar = read_string(file_outputs_from_input_t01.bar)
    String baz = read_string(file_outputs_from_input_t01.baz)
  }
}
