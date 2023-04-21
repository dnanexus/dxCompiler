
task foo {
  command {
    echo "foo"
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow wf {
  scatter (i in range(2)) {
    call foo
  }
}
