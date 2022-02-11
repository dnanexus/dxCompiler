version 1.0

# A WDL workflow with expressions, but without branches and loops.

workflow wf_linear {
  input {
    Int x = 10
  }
  call python_command {
    input: num=x
  }
  String out = python_command.result + "Did it?"
  call mul { input: a=out, b="It did NOT" }

  output {
    Int result = mul.result
  }
}


task python_command {
  input {
    Int num
  }
  command <<<
python <<CODE
mock_array = list(range(9))
failure = mock_array["${num}"]
CODE
  >>>
  output {
    String result = "This should fail with IndexError"
  }
}

task mul {
  input {
    String a
    String b
  }

  command {
    echo "I'm in the next task"
  }
  output {
    Int result = "${a}_${b}"
  }
}