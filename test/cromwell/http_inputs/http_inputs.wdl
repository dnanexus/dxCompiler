version 1.0


workflow http_inputs {
  call sum
}


task sum {

  input {
    File jamie
  }

  command {
    /usr/bin/md5sum ${jamie} | cut -c1-32
  }

  output {
    String sum = read_string(stdout())
  }

  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}
