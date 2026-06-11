version 1.0

task hello {
  input {
    String name
  }
  command <<<
    echo "Hello, ~{name}"
  >>>
  output {
    String greeting = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
