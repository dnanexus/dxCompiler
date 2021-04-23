version 1.0

workflow nested_scatter2 {
  input {
    Array[File] samples
    File reference
  }

  scatter (sample in samples) {
    call foo {
      input: sample = sample
    }
  }

  scatter (sample2 in foo.result) {
    call bar {
      input: sample = sample2, reference = reference
    }
  }

  output {
    Array[String] result = bar.s
  }
}

task foo {
  input {
    File sample
  }

  command <<<
    cat ~{sample}
  >>>

  output {
    File result = stdout()
  }
}

task bar {
  input {
    File sample
    File reference
  }

  command <<<
  cat ~{sample}
  cat ~{reference}
  >>>

  output {
    String s = read_string(stdout())
  }
}