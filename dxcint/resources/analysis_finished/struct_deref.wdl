version 1.0

struct SampleStruct {
  String sample_name
  Int id
}

workflow struct_deref {
  input {
    SampleStruct sampleStruct
    Object sampleObject
  }

  String sample_name = sampleStruct.sample_name

  call struct_deref_test {
    input:
      sample_name = sample_name,
      id = sampleStruct.id
  }

  String sample_name2 = sampleObject.sample_name

  call struct_deref_test as test2 {
    input:
      sample_name = sample_name2,
      id = sampleObject.id
  }

  output {
    String out = struct_deref_test.out
    String out2 = test2.out
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
