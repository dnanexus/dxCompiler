version 1.0

struct Result {
  String s
}

workflow nested_wf {
  input {
    String x
  }

  output {
    Result result = object {
      s: x
    }
  }
}