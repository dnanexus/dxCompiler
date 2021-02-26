version 1.0

workflow foo {
  input {
    Array[File]+ inp
  }
  output {
    Array[File]+ out = inp
  }
}