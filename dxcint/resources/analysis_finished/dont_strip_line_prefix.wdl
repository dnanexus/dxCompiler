task prefix {
  command {
    echo hello \
    | wc -l
  }
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow wf_prefix {
  call prefix
  output {
     String prefix_out = prefix.out
  }
}
