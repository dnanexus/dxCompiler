version 1.0

task yo {
  input {
    String salutation
  }
  command {
    echo 'sup ~{salutation}?'
  }

  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }

  output {
    String out = read_string(stdout())
  }
}

workflow call_cache_hit_prefixes {
  call yo

  output {
    String sup = yo.out
  }
}
