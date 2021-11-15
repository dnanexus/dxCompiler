task increment {
  Int i
  command {
    echo $(( ${i} + 1 ))
  }
  output {
    Int j = read_int(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow subwf_wt {
  Array[Int] is
  scatter (i in is) {
    call increment { input: i = i }
  }
  output {
    Array[Int] js = increment.j
  }
}
