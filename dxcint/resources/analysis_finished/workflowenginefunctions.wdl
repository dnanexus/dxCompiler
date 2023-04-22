task printArray {
  Array[Int] array

  command {
    echo "${sep='"; echo "' array}"
  }
  output {
    String o = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow workflowEngineFunctions {
  Array[Int] some_ints = range(10)
  call printArray {input: array=some_ints}
}
