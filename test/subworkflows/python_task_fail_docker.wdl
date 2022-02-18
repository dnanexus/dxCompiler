version 1.0

# A WDL workflow with a subprocess that fails with IndexError

workflow python_task_fail_docker {
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
  runtime {
    cpu: 1
    memory: "3 GB"
    docker: "python@sha256:8d125004fdd81d53a8da1cbf7d830875bbe64eb01e602b6abe9bdd5b1faffa75"
  }
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