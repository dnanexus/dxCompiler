version 1.0

workflow optional_output {
  input {
    Boolean a = true
  }
  if (a) {
    call test
    File sample = test.output_sample[0]
  }
  output {
    File? output_sample = sample
  }
}

task test {
  command {
    touch sample
  }
  output {
    Array[File] output_sample = ["sample"]
  }
}