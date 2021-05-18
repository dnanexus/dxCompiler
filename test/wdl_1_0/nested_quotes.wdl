version 1.0

task nested_quotes {
  input {
    Boolean b
    String s
  }

  command <<<
    echo ~{true='hello "~{s}"' false='goodbye "~{s}"' b}
  >>>

  output {
    String out = read_string(stdout())
  }
}