task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  runtime {
      docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    }
  output {
    String salutation = read_string(stdout())
  }
}

workflow missing_input_failure {
  File wf_hello_input
  
  call hello {input: addressee = read_string(wf_hello_input) }
  
  output {
    String salutation = hello.salutation
  }
}
