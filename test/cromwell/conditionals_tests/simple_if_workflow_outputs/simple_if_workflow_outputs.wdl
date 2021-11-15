task runMe {
  command {
    echo "done"
  }
  output {
    String s = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow simple_if {
  if (true) {
    call runMe as runMeTrue
  }

  if (false) {
    call runMe as runMeFalse
  }

  output {
    String? rmt = runMeTrue.s
    String? rmf = runMeFalse.s
  }
}
