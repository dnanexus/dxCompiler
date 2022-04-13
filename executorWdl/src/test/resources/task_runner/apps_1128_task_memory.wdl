version 1.1


task mem_int {
   input {
    Int num = 1
  }

  command <<<
    echo "a" > input_"${num}".txt
  >>>

  output {
    String result = read_string(stdout())
  }

  runtime {
    memory: "30 GB"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}
