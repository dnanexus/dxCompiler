version 1.0

struct SampleStruct {
  String sample_name
}

workflow struct_deref {
  input {
    SampleStruct sampleStruct
  }

  call test {
    input:
      sample_name = sampleStruct.sample_name
  }

  output {
    String out = test.out
  }
}

task test {
  input {
    String sample_name
  }

  command <<<
    echo ~{sample_name}
  >>>

  output {
    String out = read_string(stdout())
  }
}
