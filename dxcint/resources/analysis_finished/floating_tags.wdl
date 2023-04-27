task echo {
  command {
    echo "Peter Piper picked a peck of pickled peppers"
  }
  output {
    File out = stdout()
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

task find {
  String match = "r"
  File in_file
  command {
    grep '${match}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

workflow floating_tags {
  call echo
  call find { input: in_file = echo.out }
  call echo as echoAgain
  call find as findAgain { input: in_file = echo.out }
  output {
    Int count = find.count
    Int count2 = findAgain.count
  }
}
