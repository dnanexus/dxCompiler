task blah {
  command <<<
    python <<CODE
    def a():
      return "a"
    def b():
      return "b"
    print('{}{}'.format(a(),b()))
    CODE
  >>>

  output {
    String ab = read_string(stdout())
  }

  runtime {
    docker: "dx://file-GJ7q5KQ0yzZv7GQzJjX2x9Pb"
  }
}

workflow wf {
  call blah
}
