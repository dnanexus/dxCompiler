version 1.0

# A WDL workflow with a subprocess that fails with IndexError

workflow python_task_fail {
  call python_command
  String out = python_command.result + "Did it?"
  call mul { input: a=out, b="It did NOT" }

  output {
    String result = mul.result
  }
}

task python_command {
  command <<<
python3 <<CODE
mock_array = list(range(9))
failure = mock_array[10]
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
    String result = "${a}_${b}"
  }
}