task prepare {
 command <<<
    python3 -c "print('one\ntwo\nthree\nfour')"
  >>>
  output {
    Array[String] array = read_lines(stdout())
  }
  runtime {
    docker: "dx://file-GJ7q5KQ0yzZv7GQzJjX2x9Pb"
  }
}

task analysis {
  String str
  command <<<
    python3 -c "print('_${str}_')" > a.txt
  >>>
  output {
    File out = "a.txt"
  }
  runtime {
    docker: "dx://file-GJ7q5KQ0yzZv7GQzJjX2x9Pb"
  }
}

task gather {
  Array[File] array
  command <<<
    cat ${sep=' ' array}
  >>>
  output {
    String str = read_string(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow scattergather {
  call prepare
  scatter (x in prepare.array) {
    call analysis {input: str=x}
  }
  call gather {input: array=analysis.out}
  output {
    String gather_out = gather.str
  }
}
