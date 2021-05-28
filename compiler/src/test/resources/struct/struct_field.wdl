version 1.0

struct SampleObject {
  String s3_bucket
  Int id
}

workflow Foo {
  input {
    SampleObject sampleObject
  }

  call Baz {
    input:
      a = sampleObject.s3_bucket
  }
}

task Baz {
  input {
    String a
  }

  command <<<
    echo ~{a}
  >>>
}