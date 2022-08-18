version 1.0

workflow mock_1 {
    input {
        String in_1
    }
    call mock_1_task_1 {
        input: t1_inp = in_1
    }
    output {
        String out = mock_1_task_1.done
    }
}

task mock_1_task_1 {
  input {
    String t1_inp
  }

  command {
    echo ~{t1_inp}
  }

  output {
    String done = t1_inp
  }
}