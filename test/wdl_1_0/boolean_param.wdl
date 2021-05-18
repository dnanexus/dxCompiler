version 1.0

workflow boolean_param {
  input {
    Boolean b = true
  }

  if (b) {
    String s = "hello"
  }

  call boolean_task {
    input: b = b
  }

  output {
    String? sout = s
    String sout2 = boolean_task.s
  }
}

task boolean_task {
  input {
    Boolean b = true
  }

  command <<<
    echo ~{b}
  >>>

  output {
    String s = read_string(stdout())
  }
}