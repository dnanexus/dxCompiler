version 1.0

struct Result {
  String s
}

workflow wf {
  input {
    String x
  }

  output {
    Result result = object {
      s: x
    }
  }
}