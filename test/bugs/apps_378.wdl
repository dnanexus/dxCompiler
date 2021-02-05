version 1.0

workflow apps_378 {
  input {
    Array[String] sample_names
    Array[File] fastq1
    Array[File] fastq2
    File sample_types
  }

  Map[String, String] sample_type_map = read_map(sample_types)
  # Each row in the manifest is (sample_name, (fq1, fq2))
  # For now we're assuming it's always paired-end data
  Array[Pair[String, Pair[File, File]]] manifest = zip(sample_names, zip(fastq1, fastq2))

  scatter (row in manifest) {
    call get_sample {
      input:
        one_row = row
    }
    String sample_type = sample_type_map[get_sample.sample_name]
  }

  output {
    Array[String] sample_names_out = get_sample.sample_name
    Array[String] sample_types_out = sample_type
  }
}

task get_sample {
  input {
    Pair[String, Pair[File, File]] one_row
  }

  command {
    echo "Hello"
  }

  output {
    String sample_name = one_row.left
  }
}
