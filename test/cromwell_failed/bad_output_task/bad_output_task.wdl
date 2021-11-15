task bad {
  command {
    echo "hello" > a
  }

  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }

  output {
    # Oops! we made a spelling mistake in our WDL!
    File a = "b"
  }
}

workflow badExample {
  call bad
}
