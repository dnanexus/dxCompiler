
workflow subhello {
  Array[String] greeting_pieces

  call hello {
    input: inputs = greeting_pieces
  }

  Int salutation_length = length(hello.out)

  output {
    Int sal_len = salutation_length
  }
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
