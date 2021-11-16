import "sub_interactions_import.wdl" as counter

task hello {
  String addressee
  
  command {
    echo "Hello ${addressee}!" > hello
    wc -w < hello > count
  }
  runtime {
      docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
  output {
    String salutation = read_string("hello")
    Int count = read_int("count")
  }
}

workflow sub_workflow_interactions {
  call hello { input: addressee = "Sub Workflow World" }
  call counter.countEvens { input: max = hello.count } # Sub workflow depends on previous task call
  call hello as secondHello { input: addressee = countEvens.someStringOutput } # Task call depends on previous sub workflow call
  
  output {
    String hello_salutation = hello.salutation
    Int hello_count = hello.count

    String out2 = secondHello.salutation
    Array[Int] out4 = read_lines(countEvens.evenFile)
  }
}
