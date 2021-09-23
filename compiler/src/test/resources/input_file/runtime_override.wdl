version 1.0

task runtime_override {
  input {
    String s
  }

  command <<<
  echo ~{s}
  >>>

  runtime {
    docker: "ubuntu:latest"
    memory: "10M"
    cpu: 4
  }
}