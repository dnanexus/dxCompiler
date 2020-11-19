version 1.0

workflow archive_test {
  input {
    # workflow input is array:file
    Array[File] fa
    # workflow input is array:string
    Array[String] sa
  }

  call foo {
    input:
      # applet input is array:file
      # stage input is link from workflow input fa
      fa = fa
  }

  # applet input is array:file
  # stage input is link from workflow input fa
  scatter (p in zip(sa, fa)) {
    File g = fa
    Pair[String, File] q = p
  }

  call foo { input: fa = g }

  call bar { input: q = as_map(q) }

  output {
    Array[File] out = bar.out
  }
}

task foo {
  input {
    # task input is array:file
    Array[File] fa
  }
  command {}
  output {
    # task output: linked from task input, or is it localized and converted to an archive?
    Array[File] out = fa
  }
}

task bar {
  input {
    Map[String, File] q
  }
  command {}
  output {
    Map[String, File] out = q
  }
}