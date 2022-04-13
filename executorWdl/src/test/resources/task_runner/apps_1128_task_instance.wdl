version 1.1


task instance_type {
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
    dx_instance_type: "mem1_ssd1_v2_x16"
    dx_app: object {
      type: "app",
      id: "app-BZ9ZQzQ02VP5gpkP54b96pYY",
      name: "file_concatenator/2.0.0"
    }
  }
}
