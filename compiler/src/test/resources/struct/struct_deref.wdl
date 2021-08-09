version 1.0

struct SampleStruct {
  String sample_name
  Int id
}

workflow struct_deref {
  input {
    SampleStruct sampleStruct
  }

  String sample_name = sampleStruct.sample_name

  call struct_deref_test {
    input:
      sample_name = sample_name,
      id = sampleStruct.id
  }

  output {
    String out = struct_deref_test.out
  }
}

task struct_deref_test {
  input {
    String sample_name
    Int id
  }

  command <<<
    echo "~{sample_name} ~{id}"
  >>>

  output {
    String out = read_string(stdout())
  }
}
