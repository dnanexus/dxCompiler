version 1.1

task upload_task {
  command <<<
  set -euxo pipefail
  echo "A" > output.txt
  >>>

  output {
    File output_f = "output.txt"
  }
}

workflow upload_wf {
  call upload_task

  output {
    File output_f = upload_task.output_f
  }
}
