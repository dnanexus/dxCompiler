task make_a_file {
  command {
    echo "first line\nsecond line\nthird line" > canned
  }
  output {
    File a_file = "canned"
  }
  runtime {docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"}
}

task cat {
  File file
  String? flags

  command {
    cat ${flags} ${file}
  }

  runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }

  output {
    File procs = stdout()
  }
}

task cgrep {
  String str_decl
  String pattern
  File in_file

  command {
    grep '${pattern}' ${in_file} | wc -l
  }

  runtime { docker: "dx://file-G66qpGj0yzZq02K9313pJg5G" }

  output {
    Int count = read_int(stdout())
    String str = str_decl
  }
}

workflow two_step {
  String flags_suffix
  String flags = "-" + flags_suffix
  String static_string = "foobarbaz"

  call make_a_file

  call cat {
    input: flags = flags, file = make_a_file.a_file
  }

  call cgrep {
    input: in_file = cat.procs
  }
}
