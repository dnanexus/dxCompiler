version 1.0

workflow view_and_count {
  input {
    File bam
  }

  call slice_bam {
    input: bam = bam
  }

  scatter (slice in slice_bam.slices) {
    call count_bam {
      input: bam = slice
    }
  }

  output {
    Array[Int] counts = count_bam.count
  }
}

task slice_bam {
  input {
    File bam
  }

  command <<<
    split -l 1 -a 1 ~{bam} bam
  >>>

  output {
    Array[File] slices = glob("bam*")
  }
}

task count_bam {
  input {
    File bam
  }

  command <<<
    wc -l ~{bam} | cut -d" " -f 1
  >>>

  output {
    Int count = read_int(stdout())
  }
}
