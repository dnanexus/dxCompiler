version 1.0

task hello {
  input {
    Array[String] person
  }
  command {
    echo "hello ~{sep="," person}"
  }
  output {
    String greeting = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

task bonjour {
  input {
    Array[String]+ person
  }
  command {
    echo "hello ~{sep="," person}"
  }
  output {
    String greeting = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow postfix_quantifiers {
  input {
    Array[String] three = ["alice", "bob", "charles"]
    Array[String] one = ["alice"]
    Array[String] none = []
  }
  call hello { input: person = three }
  call hello as hello2 { input: person = one }
  call hello as hello3 { input: person = none }
  call bonjour { input: person = three }
  call bonjour as bonjour2 { input: person = one }
}
