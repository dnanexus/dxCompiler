
workflow subhello {
  Array[String] greeting_pieces

  call hello {
    input: inputs = greeting_pieces
  }

  String salutation_length = length(hello.out)

  # Neither the call output nor the declaration should be considered outputs to the subworkflow
}

task hello {
  Array[String] inputs

  command {}

  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
  output {
    Array[String] out = inputs
  }
}
