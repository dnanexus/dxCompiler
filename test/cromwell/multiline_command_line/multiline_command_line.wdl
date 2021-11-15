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
    docker: "dx://file-G66qz3Q0yzZfy6pg5q3yK3Kz"
  }
}

workflow wf {
  call blah
}
