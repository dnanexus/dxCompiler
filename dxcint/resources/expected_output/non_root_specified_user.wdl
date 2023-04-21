task nobody {
  command {
    whoami
  }

  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
    docker_user: "nobody"
  }

  output {
    String user = read_string(stdout())
  }
}

workflow woot {
  call nobody

  output {
    String nobodyUser = nobody.user
  }
}
