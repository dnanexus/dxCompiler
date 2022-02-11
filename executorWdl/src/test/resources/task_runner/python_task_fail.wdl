task python_command {
    Int num = 10
  command <<<
python <<CODE
mock_array = list(range(9))
failure = mock_array["${num}"]
CODE
  >>>
  output {
    String result = "This should fail with IndexError"
  }
}
