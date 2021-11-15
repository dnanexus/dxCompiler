task echo {
  String greeting
  String out
  Int one = 1

  command {
    echo "${greeting}" > ${out}.txt
    echo "${ one + 1 }" > ${one+1}.txt
  }
  output {
    String a_string = read_string("${out}.txt")
    Int two = read_int("${ one + 1 }.txt")
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow echo_wf {
  call echo
}

