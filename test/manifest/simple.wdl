version 1.0

workflow simple {
  input {
    String s
    File f
  }

  call task1 {
    input: s = s, f = f, n = 2
  }

  call task1 as task2 {
    input: s = task1.sout, f = task1.fout, n = 1
  }

  output {
    String sout = task2.sout
    File fout = task2.fout
  }
}

task task1 {
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