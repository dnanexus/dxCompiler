
workflow subhello {
  Array[String] greeting_pieces

  call hello {
    input: inputs = greeting_pieces
  }

  # Confirm referencing call outputs from subworkflows does not fail validation.
  Int salutation_length = length(hello.out)

  output {
    Array[String] hello_out = hello.out
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
