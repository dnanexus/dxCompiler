task runMe {
  command {
    echo "done"
  }
  output {
    String s = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow simple_if {
  if (true) {
    call runMe as runMeTruee
  }

  if (false) {
    call runMe as runMeFalse
  }

  output {
    String? true_res = runMeTruee.s
    String? false_res = runMeFalse.s
  }
}
