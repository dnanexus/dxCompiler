version 1.0

workflow localizer_workflow {
  input {
    File input_file
  }

  call localizer_task {
    input:
      input_file = input_file
  }

  output {
    String o = localizer_task.out
  }
}

task localizer_task {
  input {
    File input_file
  }

  command <<<
    pwd && ls -lah .
  >>>

  runtime {
  	docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }

  output {
    String out = read_string(stdout())
  }
}
