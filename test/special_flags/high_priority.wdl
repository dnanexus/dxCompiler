version 1.0

task echo {
  input {
    String s
  }
  command <<<
  echo ~{s}
  >>>
  output {
    String out = read_string(stdout())
  }
}

workflow high_priority {
  input {
    String s
  }
  call echo as echo_priority {
    input: s = s
  }
  output {
    String out = echo_priority.out
  }
}