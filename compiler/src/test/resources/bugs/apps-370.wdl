version 1.0

workflow test {
  input {
    ## Manifest: string: sample name; [File, File]: paired fastq files
    Array[Pair[String, Pair[File, File]]] manifest

    ## Each row in the tsv sample_types file is (sample_name\tsample_type)
    File sample_types
    Map[String, String] sample_type_map = read_map(sample_types)
  }

  scatter (item in manifest) {
    call testtask {
      input:
        fastq1 = item.right.left,
        fastq2 = item.right.right,
    }

    String sample_type = sample_type_map[item.left]
  }

  output {
    Array[String] out = testtask.out
    Array[String] sample_types = sample_type
  }
}

task testtask {
  input {
    File fastq1
    File fastq2
  }

  command {
    echo "hello bam"
  }

  output {
    String out = "hello"
  }
}
