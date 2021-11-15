task a {
  command {
    echo "12"
    >&2 echo "200"
  }
  output {
    File out = stdout()
    File err = stderr()
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

task b {
  File in_file
  command {
    cat ${in_file}
  }
  output {
    Int out = read_int(stdout())
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

workflow stdout_stderr_passing {
  call a
  call b {input: in_file=a.out}
  call b as b_prime {input: in_file=a.err}
  output {
    Int b_prime_out = b_prime.out
  }
}
