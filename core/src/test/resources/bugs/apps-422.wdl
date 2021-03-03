version 1.0

struct Foo {
  String s
}

workflow struct_test {
  input {
    String? wf_input
  }

  if (defined(wf_input)) {
    Foo a = {"s": "no_input"}
  }
  if (defined(wf_input)) {
    Foo b = {"s": "one_input"}
  }

  Foo c = select_first([a,b])

  call test {
    input:
      answer = c.s
  }

  output {
    String out = test.out
  }
}

task test {
  input {
    String answer
  }

  command <<<
    echo ~{answer}
  >>>

  output {
    String out = read_string(stdout())
  }
}
