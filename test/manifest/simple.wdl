version 1.0

workflow simple {
  input {
    String s
    File f
  }

  call echocat {
    input: s = s, f = f, n = 2
  }

  call echocat as echocat2 {
    input: s = echocat.sout, f = echocat.fout, n = 1
  }

  call merge {
    input:
      strings = [echocat.sout, echocat2.sout],
      files = [echocat.fout, echocat2.fout]
  }

  output {
    String sout = merge.sout
    File fout = merge.fout
  }
}

task echocat {
  input {
    String s
    File f
    Int n
  }

  command <<<
  echo "hello ~{s}"
  cat ~{f} | head -~{n} > out
  >>>

  output {
    String sout = read_string(stdout())
    File fout = "out.${n}"
  }
}

task merge {
  input {
    Array[String] strings
    Array[File] files
  }

  command <<<
    cat ~{sep=" " files} > out_merge.txt
  >>>

  output {
    String sout = "${sep=" " strings}"
    File fout = "out_merge.txt"
  }
}
