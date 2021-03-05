version 1.1

task add2 {
  input {
    Int a
    Int b
  }

  command <<<
    echo $((~{a} + ~{b}))
  >>>

  runtime {
    container: "ubuntu:latest"
  }

  meta {
    summary: "Adds two int together"
  }

  output {
    Int result = read_int(stdout())
  }
}
