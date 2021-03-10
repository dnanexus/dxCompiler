version 1.0

task diff_stream_and_download {
  input {
    File a
    File b
  }
  parameter_meta {
    a : "stream"
  }
  runtime {
    docker: "ubuntu:16.04"
  }
  command {
    diff ${a} ${b} | wc -l
  }
  output {
    Int result = read_int(stdout())
  }
}