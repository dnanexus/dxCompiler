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

workflow wf_hello {
  String wf_hello_input = "world"
  
  call hello { input: addressee = wf_hello_input }
  
  output {
    String salutation = hello.salutation
  }
}
