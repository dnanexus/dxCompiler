version 1.1

import "structs.wdl"

workflow test_workflow {
  input {
    MyStructDB database
    String prefix
  }

  call test_task2 {
    input:
      prefix = prefix
  }
  
  output {
    File final_output = test_task2.task_output
  }
}

task test_task2 {
  input {
    String prefix
  }

  command <<<
     echo ~{prefix} > output.txt
  >>>

  runtime {
    docker: "debian:stretch-slim"
  }

  output {
    File task_output = "output.txt"
  }
}
