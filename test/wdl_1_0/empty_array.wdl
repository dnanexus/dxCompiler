version 1.0

workflow empty_array {
  input {
    Array[String] names
  }

  scatter (name in names) {
    String greet = "hello ${name}"
  }

  output {
    Array[String] greetings = greet
    Int num_names = length(names)
  }
}
