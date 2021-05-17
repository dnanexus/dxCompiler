version 1.1

task task_1 {
  command <<<
  set -euxo pipefail
  echo "A" > output.txt
  >>>

  output {
    File output_f = "output.txt"
  }
}

workflow workflow_1 {
  call task_1

  output {
    File output_f = task_1.output_f
  }
}
