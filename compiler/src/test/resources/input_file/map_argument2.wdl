version 1.0

workflow map_argument2 {
  input {
    Map[String, File] files
  }
  output {
    Map[String, File] result = files
  }
}