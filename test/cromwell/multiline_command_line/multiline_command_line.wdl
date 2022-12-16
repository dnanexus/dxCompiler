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
    docker: "dx://file-GJ941b80yzZvGbK68zxQzB0B"
  }
}

workflow multiline_command_line {
  call blah
}
